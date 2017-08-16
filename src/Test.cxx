// Test the NMRMath library 

#include <ctime> 
#include <cstdlib> 
#include <iostream> 
#include <fstream> 
#include <vector>
#include <iomanip>      // std::setprecision

#include "NMRMath.hh"

template < typename T1, typename T2 > int CalculateFrequencyZC(int inputUnits,double SampleFreq,std::vector<T1> tCross,T2 &freq_full_range); 
template < typename T1, typename T2 > int PrintToFile(std::string outpath,std::vector<T1> x,std::vector<T2> y); 

int main(int argc,char **argv){

   if(argc<2){
      std::cout << "Usage: ./bin/Test <frequency (kHz)>" << std::endl;
      return 1;
   }

   double inFreq = atof(argv[1])*1E+3;

   std::vector<unsigned long> time; 
   std::vector<double> voltage; 

   double SAMPLE_FREQ = 10E+6;
   double omega       = 2.*acos(-1)*inFreq;
   double ampl        = 1.0; 

   unsigned long sample_num; 
   double arg_v,arg_t;
 
   const int NPTS = 1E+6;
   for (int i=0;i<NPTS;i++) {
      sample_num = i;     
      arg_t      = (double)(i)/SAMPLE_FREQ;
      arg_v      = ampl*sin(omega*arg_t);
      time.push_back(sample_num);  
      voltage.push_back(arg_v); 
   }

   std::string prefix  = "/Users/dflay/work/E989/NMRMath/ana/"; 
   std::string outpath = prefix + "test-data.csv"; 
   PrintToFile<unsigned long,double>(outpath,time,voltage);

   int verbosity = 0;
   int method    = NMRMath::kLinearInterpolation;

   bool UseTimeRange  = true; 
   unsigned long tMin = 0;
   unsigned long tMax = 0.75*NPTS;

   std::vector<unsigned long> tCross;
   std::vector<double> vCross; 

   const int NC = 1E+5; 
   tCross.reserve(NC);
   vCross.reserve(NC);
   for(int i=0;i<NC;i++){
      tCross.push_back(0); 
      vCross.push_back(-1); 
   }

   int start_time = clock(); 

   int numCrossings = NMRMath::CountZeroCrossings(verbosity,method,SAMPLE_FREQ,inFreq,UseTimeRange,tMin,tMax,time,voltage,tCross,vCross);

   std::cout << "Found " << numCrossings << " zero crossings " << std::endl;

   // clean up invalid points from initialization 
   for(int i=NC-1;i>=numCrossings;i--){
      tCross.pop_back();
      vCross.pop_back();
   }

   std::string outpath_zc = prefix + "zc-calc.csv"; 
   PrintToFile<unsigned long,double>(outpath_zc,tCross,vCross); 
 
   double outFreq=0; 
   int rc = CalculateFrequencyZC<unsigned long,double>(NMRMath::kNumSamples,SAMPLE_FREQ,tCross,outFreq);

   int end_time = clock();

   double duration = (end_time-start_time)/( (double)CLOCKS_PER_SEC );  

   char inf[20],outf[20],diff[20],dur[20];

   double df = (outFreq-inFreq)/0.06179;
 
   sprintf(inf ,"%.5lf Hz",inFreq );
   sprintf(outf,"%.5lf Hz",outFreq);
   sprintf(diff,"%.5lf ppb",df);
   sprintf(dur ,"%.3lf sec",duration); 

   std::cout << "Input frequency:  " << inf  << std::endl; 
   std::cout << "Output frequency: " << outf << std::endl; 
   std::cout << "Difference:       " << diff << std::endl; 
   std::cout << "Elapsed time:     " << dur  << std::endl;  
 
   return 0;
}
//______________________________________________________________________________
template < typename T1, typename T2 > 
int CalculateFrequencyZC(int inputUnits,double SampleFreq,std::vector<T1> tCross,T2 &freq_full_range){

   int rc = 0;
   int NumCrossings = tCross.size();
   std::cout << "Number of crossings = " << NumCrossings << std::endl; 

   if(NumCrossings==0) std::cout << "[NMRZeroCrossing::CalculateFrequencies]: NumCrossings is zero!" << std::endl;

   T2 Zc = (T2)NumCrossings;
   T2 NC = ( Zc - (T2)1 )/( (T2)2 );         // number of cycles

   if(NumCrossings<5){
      std::cout << "[NMRZeroCrossings::CalculateFrequencies]: ERROR!  Number of zero crossings is less than 5! " << std::endl;
      freq_full_range = 0;
      return 1;
   }

   // find frequency using all zero crossings 
   // find the index of last crossing  
   int LastCrossingIndex  = NumCrossings - 1;  // the -1 is to account for arrays starting at 0 and not 1

   // compute the time 
   T2 dt  = (T2)( tCross[LastCrossingIndex] - tCross[0] );
   // compute the frequency  
   freq_full_range = NC/dt;

   if(inputUnits==NMRMath::kNumSamples) {
      freq_full_range *= SampleFreq;
   }

   return rc;
}
//______________________________________________________________________________
template < typename T1, typename T2 > 
int PrintToFile(std::string outpath,std::vector<T1> x,std::vector<T2> y){

   const int N = x.size();

   std::ofstream outfile;
   outfile.open(outpath); 
   if ( outfile.fail() ) {
      std::cout << "Cannot write the data to the file: " << outpath << std::endl;
      return 1; 
   } else { 
      for (int i=0;i<N;i++) outfile << x[i] << "," << y[i] << std::endl;
      outfile.close();
      std::cout << "The data has been written to the file: " << outpath << std::endl;
   }

   return 0;  
}
