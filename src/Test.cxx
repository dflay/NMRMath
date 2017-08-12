// Test the NMRMath library 

#include <cstdlib> 
#include <iostream> 
#include <vector> 

#include "NMRMath.hh"

int main(){

   std::cout << "Vector contains: " << std::endl;
   std::vector<double> A;
   A.push_back(6);
   A.push_back(7);
   A.push_back(9);
   const int N = A.size();
   for(int i=0;i<N;i++) std::cout << A[i] << std::endl;

   double mean = NMRMath::GetMean<double>(A);
   std::cout << "mu = " << mean << std::endl;

   double stdev = NMRMath::GetStandardDeviation<double>(A);
   std::cout << "sigma = " << stdev << std::endl;

   return 0;
}
