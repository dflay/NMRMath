#ifndef NMRMATH_HH
#define NMRMATH_HH

#include <cstdlib>
#include <iostream>
#include <vector>
#include <math.h>

namespace NMRMath { 

   enum zcType{
      kMidpoint            = 1,
      kLinearInterpolation = 2,
      kLeastSquares        = 3
   }; 

   //______________________________________________________________________________
   template < class T >
      T GetMean(std::vector<T> x){
         int N = x.size();
         T sum=0;
         for(int i=0;i<N;i++) sum += x[i];
         T mean = sum/( (T)N );
         return mean;
      }
   //______________________________________________________________________________
   template < class T >
      T GetRMS(std::vector<T> x){
         T sum_sq = 0;
         const int N = x.size();
         for(int i=0;i<N;i++){
            sum_sq += pow(x[i],2.);
         }
         T arg = sum_sq/( (T)N );
         T rms = sqrt(arg);
         return rms;
      }
   //______________________________________________________________________________
   template < class T >
      T GetVariance(std::vector<T> x){
         int N = x.size();
         T mean = GetMean<T>(x);
         T sum=0;
         for(int i=0;i<N;i++){
            sum += pow(x[i]-mean,2);
         }
         T var   = sum/( (T)N );
         return var;
      }
   //______________________________________________________________________________
   template < class T >
      T GetStandardDeviation(std::vector<T> x){
         T var   = GetVariance<T>(x);
         T stdev = sqrt(var);
         return stdev;
      }
   //______________________________________________________________________________
   template < class T >
      T GetStandardErrorOfTheMean(std::vector<T> x){
         const int N = x.size();
         T sd   = GetStandardDeviation<T>(x);
         T sdom = sd/sqrt( (T)N );
         return sdom;
      }
   //______________________________________________________________________________
   template < class T >
      T GetCovariance(std::vector<T> x,std::vector<T> y){
         T mean_x = GetMean<T>(x);
         T mean_y = GetMean<T>(y);
         T sum=0,diff_x=0,diff_y=0;
         const int N = x.size();
         for (int i=0;i<N;i++) {
            diff_x = x[i]-mean_x;
            diff_y = y[i]-mean_y;
            sum   += diff_x*diff_y;
         }
         T cov = sum/( (T)N );
         return cov;
      }
   //______________________________________________________________________________
   template < class T >
      int LeastSquaresFitting(std::vector<T> x,std::vector<T> y,T &a,T &b,T &r){
         // linear regression to find slope b and y-intercept a of 
         // f(x) = a + bx 

         int rc=0;
         double num=0,rsq=0;

         int N    = x.size();
         T mean_x = GetMean<T>(x);
         T mean_y = GetMean<T>(y);
         T var_x  = GetVariance<T>(x);
         T var_y  = GetVariance<T>(y);
         T cov_xy = GetCovariance<T>(x,y);

         T ss_xx = ( (T)N )*var_x;
         T ss_yy = ( (T)N )*var_y;
         T ss_xy = ( (T)N )*cov_xy;

         T den = ss_xx*ss_yy;
         if(den==0){
            // singular matrix. can't solve the problem.
            a   = 0;
            b   = 0;
            r   = 0;
            rc  = 1;
         }else{
            b   = cov_xy/var_x;
            a   = mean_y - b*mean_x;
            num = ss_xy*ss_xy;
            rsq = num/den;
            r   = sqrt(rsq);
         }

         return rc;
      }

} // ::NMRMath 

#endif 
