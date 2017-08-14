// Test the NMRMath library 

#include <cstdlib> 
#include <iostream> 
#include <vector> 

#include "NMRMath.hh"

template < typename T1, typename T2 > int CalculateFrequencyZC(std::vector<T1> tCross,T2 &freq_full_range); 

int main(){

   std::vector<double> A;
   A.push_back(6);
   A.push_back(7);
   A.push_back(9);
   const int N = A.size();

   // std::cout << "Vector contains: " << std::endl;
   // for(int i=0;i<N;i++) std::cout << A[i] << std::endl;

   double mean = NMRMath::GetMean<double>(A);
   std::cout << "mu = " << mean << std::endl;

   double stdev = NMRMath::GetStandardDeviation<double>(A);
   std::cout << "sigma = " << stdev << std::endl;

   // try interpolation 
   std::vector<unsigned long> B; 
   B.push_back(1000); 
   B.push_back(1400); 
   B.push_back(3000); 
   // for(int i=0;i<N;i++) std::cout << time[i] << std::endl;

   unsigned long myTime = B[0] + 200;

   double A_pt = NMRMath::LinearInterpolation<unsigned long,double>(myTime,B[0],A[0],B[1],A[1]);
   std::cout << myTime << " " << A_pt << std::endl;

   // now try zero crossings 
   std::vector<unsigned long> time; 
   std::vector<double> voltage; 

   double inFreq = 10E+3; 
   double SAMPLE_FREQ = 10E+6;
   double omega = 2.*acos(-1)*inFreq;
   double ampl = 1.0; 

   unsigned long sample_num; 
   double arg_v,arg_t;
 
   const int NPTS = 1E+5; 
   for (int i=0;i<NPTS;i++) {
      sample_num = i;     
      arg_t = (double)(i)/SAMPLE_FREQ;
      arg_v = ampl*sin(omega*arg_t);
      time.push_back(sample_num);  
      voltage.push_back(arg_v); 
      // std::cout << i << "  " << arg_t << "  " << arg_v << std::endl; 
   }

   int verbosity = 5;
   int method    = NMRMath::kLinearInterpolation;

   double N_exp  = SAMPLE_FREQ/inFreq;       // number of points for one period 
   int step_size = (int)( (1./16.)*N_exp );  // skip 1/16 of a period 
   int npts      = step_size/2;              // use step_size/2 

   bool UseTimeRange = false; 
   unsigned long tMin = 0;
   unsigned long tMax = 1;

   std::vector<unsigned long> tCross;
   std::vector<double> vCross; 

   int rc = NMRMath::CountZeroCrossings(verbosity,method,npts,step_size,UseTimeRange,tMin,tMax,time,voltage,tCross,vCross);

   double outFreq=0; 
   rc = CalculateFrequencyZC<unsigned long,double>(tCross,outFreq);
   outFreq *= SAMPLE_FREQ;  // because the time is sample number, we have to convert to Hz  

   std::cout << "Input frequency:  " << inFreq << std::endl; 
   std::cout << "Output frequency: " << outFreq << std::endl; 
 
   return 0;
}
//______________________________________________________________________________
template < typename T1, typename T2 > 
int CalculateFrequencyZC(std::vector<T1> tCross,T2 &freq_full_range){

   int rc = 0;
   int NumCrossings = tCross.size();

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

   return rc;
}
