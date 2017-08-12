// Test the NMRMath library 

#include <cstdlib> 
#include <iostream> 
#include <vector> 

#include "NMRMath.hh"

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
   std::vector<unsigned int> time; 
   time.push_back(1000); 
   time.push_back(1400); 
   time.push_back(3000); 
   // for(int i=0;i<N;i++) std::cout << time[i] << std::endl;

   unsigned int myTime = time[0] + 200;

   double A_pt = NMRMath::LinearInterpolation<unsigned int,double>(myTime,time[0],A[0],time[1],A[1]);
   std::cout << myTime << " " << A_pt << std::endl;
 
   return 0;
}
