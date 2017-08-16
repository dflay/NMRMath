#ifndef NMRMATH_HH
#define NMRMATH_HH

#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath> 

namespace NMRMath { 
   //______________________________________________________________________________
   enum zcType{
      kMidpoint            = 1,
      kLinearInterpolation = 2,
      kLeastSquares        = 3
   }; 
   //______________________________________________________________________________
   enum inputUnits{
      kNumSamples   = 1,
      kSeconds      = 2,
      kMilliSeconds = 3,
      kMicroSeconds = 4,
      kNanoSeconds  = 5
   }; 
   //______________________________________________________________________________
   template < typename T >
      T GetMean(std::vector<T> x){
         int N = x.size();
         T sum=0;
         for(int i=0;i<N;i++) sum += x[i];
         T mean = sum/( (T)N );
         return mean;
      }
   //______________________________________________________________________________
   template < typename T >
      T GetMean(const int N,std::vector<T> x){
         T sum=0;
         for(int i=0;i<N;i++) sum += x[i];
         T mean = sum/( (T)N );
         return mean;
      }
   //______________________________________________________________________________
   template < typename T >
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
   template < typename T >
      T GetRMS(const int N,std::vector<T> x){
         T sum_sq = 0;
         for(int i=0;i<N;i++){
            sum_sq += pow(x[i],2.);
         }
         T arg = sum_sq/( (T)N );
         T rms = sqrt(arg);
         return rms;
      }
   //______________________________________________________________________________
   template < typename T >
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
   template < typename T >
      T GetVariance(const int N,std::vector<T> x){
         T mean = GetMean<T>(N,x);
         T sum=0;
         for(int i=0;i<N;i++){
            sum += pow(x[i]-mean,2);
         }
         T var   = sum/( (T)N );
         return var;
      }
   //______________________________________________________________________________
   template < typename T >
      T GetStandardDeviation(std::vector<T> x){
         T var   = GetVariance<T>(x);
         T stdev = sqrt(var);
         return stdev;
      }
   //______________________________________________________________________________
   template < typename T >
      T GetStandardDeviation(const int N,std::vector<T> x){
         T var   = GetVariance<T>(N,x);
         T stdev = sqrt(var);
         return stdev;
      }
   //______________________________________________________________________________
   template < typename T >
      T GetStandardErrorOfTheMean(std::vector<T> x){
         const int N = x.size();
         T sd   = GetStandardDeviation<T>(x);
         T sdom = sd/sqrt( (T)N );
         return sdom;
      }
   //______________________________________________________________________________
   template < typename T >
      T GetStandardErrorOfTheMean(const int N,std::vector<T> x){
         T sd   = GetStandardDeviation<T>(N,x);
         T sdom = sd/sqrt( (T)N );
         return sdom;
      }
   //______________________________________________________________________________
   template < typename T1, typename T2 >
      T2 GetCovariance(std::vector<T1> x,std::vector<T2> y){
         T1 mean_x = GetMean<T1>(x);
         T2 mean_y = GetMean<T2>(y);
         T2 sum=0,diff_x=0,diff_y=0;
         const int N = x.size();
         for (int i=0;i<N;i++) {
            diff_x = (T2)(x[i]-mean_x);
            diff_y = y[i]-mean_y;
            sum   += diff_x*diff_y;
         }
         T2 cov = sum/( (T2)N );
         return cov;
      }
   //______________________________________________________________________________
   template < typename T1, typename T2 >
      T2 GetCovariance(const int N,std::vector<T1> x,std::vector<T2> y){
         T1 mean_x = GetMean<T1>(N,x);
         T2 mean_y = GetMean<T2>(N,y);
         T2 sum=0,diff_x=0,diff_y=0;
         for (int i=0;i<N;i++) {
            diff_x = (T2)(x[i]-mean_x);
            diff_y = y[i]-mean_y;
            sum   += diff_x*diff_y;
         }
         T2 cov = sum/( (T2)N );
         return cov;
      }
   //______________________________________________________________________________
   template < typename T1, typename T2 >
      int LeastSquaresFitting(std::vector<T1> x,std::vector<T2> y,T2 &a,T2 &b,T2 &r){
         // linear regression to find slope b and y-intercept a of 
         // f(x) = a + bx 

         int rc=0;
         double num=0,rsq=0;

         const int N = x.size(); 
         T1 mean_x = GetMean<T1>(x);
         T2 mean_y = GetMean<T2>(y);
         T1 var_x  = GetVariance<T1>(x);
         T2 var_y  = GetVariance<T2>(y);
         T2 cov_xy = GetCovariance<T1,T2>(x,y);

         T1 ss_xx = ( (T1)N )*var_x;
         T2 ss_yy = ( (T2)N )*var_y;
         T2 ss_xy = ( (T2)N )*cov_xy;

         T2 den = (T2)(ss_xx)*ss_yy;
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
   //______________________________________________________________________________
   template < typename T1, typename T2 >
      int LeastSquaresFitting(const int N,std::vector<T1> x,std::vector<T2> y,T2 &a,T2 &b,T2 &r){
         // linear regression to find slope b and y-intercept a of 
         // f(x) = a + bx 

         int rc=0;
         double num=0,rsq=0;

         T1 mean_x = GetMean<T1>(N,x);
         T2 mean_y = GetMean<T2>(N,y);
         T1 var_x  = GetVariance<T1>(N,x);
         T2 var_y  = GetVariance<T2>(N,y);
         T2 cov_xy = GetCovariance<T1,T2>(N,x,y);

         T1 ss_xx = ( (T1)N )*var_x;
         T2 ss_yy = ( (T2)N )*var_y;
         T2 ss_xy = ( (T2)N )*cov_xy;

         T2 den = (T2)(ss_xx)*ss_yy;
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
   //______________________________________________________________________________
   template < typename TA, typename TB >
      TA LinearInterpolationForTime(TB y,TA x0,TB y0,TA x1,TB y1){
         // time can't be negative here, so we take the absolute values 
         // when calculating the slope  
	 TA b = std::abs(y-y0)/std::abs(y1-y0);
	 TA x = x0 + b*(x1-x0);
	 return x;
      }
   //______________________________________________________________________________
   template < typename TA, typename TB >
      TB LinearInterpolation(TA x,TA x0,TB y0,TA x1,TB y1){
	 TB b = (TB)(x-x0)/(TB)(x1-x0);
	 TB y = y0 + b*(y1-y0); 
	 return y;
      }
   //______________________________________________________________________________
   template <typename T1,typename T2>
   T1 GetTimeOfCrossing(int verbosity,int method,const int NN,std::vector<T1> X,std::vector<T2> Y,
                        T1 t_current,T2 v_current,
                        T1 t_next   ,T2 v_next){
      // find the time of the zero crossing 
      // t_current,v_current, etc: (time,voltage) for midpoint method 
      // X,Y = (time,voltage) to use in the fitting procedure if doing linear interpolation or least squares 

      // FIXME: We're systematically low by 1 tick mark when using T1 = unsigned int. 

      int ret_val=0;
      T1 t0=0;
      T2 v0=0,a=0,b=0,r=0;

      if (method==kMidpoint){
         // method 1: take midpoint between t_current and t_next 
         t0 = (t_current + t_next)/2.;
      } else if (method==kLinearInterpolation){
         // method 2: get time at V = 0, linear interpolation  
         t0 = LinearInterpolationForTime<T1,T2>(v0,t_current,v_current,t_next,v_next);
      } else if(method==kLeastSquares){
         // method 3: least squares fit to neighboring points
         // to find fit parameters a and b in f(x) = a + bx 
         ret_val = LeastSquaresFitting<T1,T2>(NN,X,Y,a,b,r);
         if (ret_val!=0) {
	    std::cout << "[NMRMath::GetTimeOfCrossing]: ERROR!  Least squares fit failed!" << std::endl;
	    t0 = -1;
	 } else {
	    if(b!=0){
	       t0 = (T1)( -a/b );
	    }else{
	       t0 = 0;
	    }
	 }
         // make sure t0 is bound properly 
         if(t0<X[0] || t0>X[NN-1]){
            if(verbosity>=3){
               std::cout << "[NMRMath::GetTimeOfCrossing]: ERROR!  t0 is out of bounds!  ";
               std::cout << "  Using linear interpolation instead..." << std::endl;
               std::cout << "                              t_min = "  << X[0]       << "\t" 
                                                                      << "t_max = " << X[NN-1] << "\t" 
                                                                      << "t0 = "    << t0        << std::endl;
            }
	    t0 = LinearInterpolationForTime<T1,T2>(v0,t_current,v_current,t_next,v_next);
            if(verbosity>=3) {
               std::cout << "[NMRMath::GetTimeOfCrossing]: linear interpolation: t_current = " << t_current << "\t"
                         << " t_next = " << t_next << "\t" << "t0 = " << t0 << std::endl;
            }
         }
      }else{
         // invalid method, set to -1
         t0 = -1;
      }

      if(verbosity>=3){
         if(t0<0){
            std::cout << "[NMRMath::GetTimeOfCrossing]: BAD CROSSING TIME!" << std::endl;
            std::cout << "                              t0        = " << t0        << std::endl;
            std::cout << "                              method    = " << method    << std::endl;
            std::cout << "                              t_current = " << t_current << std::endl;
            std::cout << "                              t_next    = " << t_next    << std::endl;
            if(method==kLeastSquares){
               std::cout << "                              offset = " << a << std::endl;
               std::cout << "                              slope  = " << b << std::endl;
               for(int i=0;i<NN;i++){
        	  std::cout << "                              t = " << X[i] << "\t" << "v = " << Y[i] << std::endl;
               }
            }
         }
      }

      return t0;
   }
   //______________________________________________________________________________
   template < typename T1, typename T2 > 
   int InitializeVectors(size_t N,std::vector<T1> &X,std::vector<T2> &Y){
      X.reserve(N);
      Y.reserve(N);
      for(size_t i=0;i<N;i++){
	 X.push_back(0);
	 Y.push_back(0);
      }
      return 0;
   }
   //______________________________________________________________________________
   template < typename T1, typename T2 > 
   int ClearVectors(std::vector<T1> &X,std::vector<T2> &Y){
      const int N = X.size(); 
      for(int i=0;i<N;i++){
	 X[i] = 0;
	 Y[i] = 0;
      }
      return 0;
   }
   //______________________________________________________________________________
   template < typename T1, typename T2> 
   int StoreData(int verbosity,int i,int NPTS,
                 std::vector<T1> time,std::vector<T2> voltage,
                 std::vector<T1> &X,std::vector<T2> &Y){

      int N     = time.size();
      int start = i - NPTS/2;
      int end   = i + NPTS/2;

      T2 v_current = voltage[i];

      if(v_current!=0){
         // do nothing
      }else{
         // voltage of the zero crossing is zero! 
         start = i - 3;
         end   = i + 3;
         if(verbosity>=3){
            std::cout << "[NMRMath::StoreData]: Voltage at zero crossing is zero!" << std::endl;
            std::cout << "                      start = " << start << std::endl;
            std::cout << "                      end   = " << end   << std::endl;
         }
      }

      // prevent unrealistic bounds: use 6 data points for the fit 
      if(start < 0){
         start = i;
         end   = i + 6;
         if(verbosity>=3){
            std::cout << "[NMRMath::StoreData]: Invalid start point!  Setting to index = " << i << std::endl;
            std::cout << "                      start = " << start << std::endl;
            std::cout << "                      end   = " << end   << std::endl;
         }
      }else if(start==0){
         start = 0;
         end   = 7;
         if(verbosity>=3){
            std::cout << "[NMRMath::StoreData]: Starting at index = " << start << std::endl;
            std::cout << "                      start = " << start << std::endl;
            std::cout << "                      end   = " << end   << std::endl;
         }
      }

      if(end > N){
	 start = N - NPTS;
	 end   = N;
      }
      int k=0;
      for(int j=start;j<end;j++){
	 X[k] = time[j]; 
	 Y[k] = voltage[j];
	 k++;
      }

      int NPTSUseable= NPTS;

      if(k!=NPTS){
	 NPTSUseable = k;
	 if(verbosity>=3){
	    std::cout << "[NMRMath::StoreData]: WARNING!  Do not have the expected number of data points! " << std::endl;
	    std::cout << "                      k    = " << k    << std::endl;
	    std::cout << "                      NPTS = " << NPTS << std::endl;
	    std::cout << "                      Using k data points in fit..." << std::endl;
	 }
      }

      return NPTSUseable;
   } 
   //______________________________________________________________________________
   template < typename T1, typename T2 >
   int CountZeroCrossings(int verbosity,int method,double SampleFreq,double ExpFreq,
                          bool UseTimeRange,T1 tMin,T1 tMax,
                          std::vector<T1> time,std::vector<T2> voltage, 
                          std::vector<T1> &tCross,std::vector<T2> &vCross){
      // count zero crossings  
      // verbosity:    how much info is printed to screen
      // method:       midpoint, linear interpolation, least squares 
      // UseTimeRange: use a specific time range specified by tMin, tMax 
      // time:         time vector 
      // voltage:      voltage vector 
      // tCross:       time of zero crossing 
      // vCross:       voltage of zero crossing (should be zero) 

      double N_exp  = SampleFreq/ExpFreq;       // number of points for one period 
      int step_size = (int)( (1./16.)*N_exp );  // skip 1/16 of a period 
      int NPTS      = step_size/2;              // use step_size/2 
      
      // NPTS:         number of points to use in zc counting (in linear fit)  
      // step_size:    number of points to to skip after a zero crossing is found     

      if(verbosity>=3) std::cout << "[NMRMath]: Counting zero crossings..." << std::endl;

      int NPTSUseable = 0;
      const int N     = time.size();

      // for counting how many zero crossings we have 
      int cntr        = 0;
      int cntr_prev   = 0;

      T1 t0=0,t_current=0,t_previous=0,t_next=0;

      T2 v0=0,target=0;
      T2 v_prod=0,delta_v=0;
      T2 v_current=0,v_previous=0,v_next=0;

      int index=0;  // vector index
      std::vector<T1> X; 
      std::vector<T2> Y; 

      InitializeVectors<T1,T2>(NPTS,X,Y); 

      int i=0;
      do { 
         t_current     = time[i];      
         t_next        = time[i+1];    
         v_current     = voltage[i];   
         v_next        = voltage[i+1]; 
         v_prod        = v_current*v_next;
         if (v_prod>target) {
            // positive number, no crossing
            // increment i by 1
            i++;
         } else if (v_prod<=target) { 
            // negative number or ZERO; we had a crossing
            if (UseTimeRange) {
               // use the fit range
               if (t_current>tMin && t_next<tMax) {
        	  cntr++;  // count the crossing 
                  delta_v = fabs(v_current-v_next); 
                  // fill vectors for fit method 
                  NPTSUseable = StoreData<T1,T2>(verbosity,i,NPTS,time,voltage,X,Y); 
                  t0 = GetTimeOfCrossing<T1,T2>(verbosity,method,NPTSUseable,X,Y,t_current,v_current,t_next,v_next); 
                  // fill vectors 
                  tCross[index] = t0; 
                  vCross[index] = v0; 
		  index++;  
               }
            } else {
               // don't use the fit range
               cntr++; // count the crossing 
               delta_v = fabs(v_current-v_next); 
               // fill vectors for fit method 
               NPTSUseable = StoreData<T1,T2>(verbosity,i,NPTS,time,voltage,X,Y); 
               t0 = GetTimeOfCrossing<T1,T2>(verbosity,method,NPTSUseable,X,Y,t_current,v_current,t_next,v_next); 
               // fill vectors 
               tCross[index] = t0; 
               vCross[index] = v0;
               index++;  
            }
            // std::cout << " CROSSING " << cntr << std::endl; 
            // check the t0 
            // std::cout << "Checking the crossing" << std::endl; 
            if(t0<0 || delta_v>0.100){   // we shouldn't see a 100 mV jump during the zero crossing 
               if(verbosity>=4){
        	  std::cout << "[NMRMath::CountZeroCrossings]: bad crossing for Zc = "
        	            << cntr << "!  Trying next crossing..." << std::endl;
        	  std::cout << "                               t0      = " << t0      << std::endl;
        	  std::cout << "                               delta_v = " << delta_v << std::endl;
               }
               cntr--;      // don't count the crossing, decrement cntr 
               i += step_size;   // move to next bin 
               // go back an entry in the dummy vectors since we have a bad crossing   
               index--; 
               // clear analysis vectors 
               // ClearVectors<T1,T2>(X,Y);
               if(verbosity>=4){
        	  std::cout << "[NMRMath::CountZeroCrossings]: zero crossing counter reset to: " << cntr << std::endl;
        	  std::cout << "[NMRMath::CountZeroCrossings]: moving to index:                " << i       << std::endl;
               }
               continue;
            }
            // passed t0 check 
            // std::cout << "CROSSING " << cntr << "  " << t0 << std::endl;
            i += step_size; // move to next bin  
            // set up for next data point 
	    // ClearVectors<T1,T2>(X,Y); 
            // index          = 0; 
            cntr_prev      = cntr; 
            t_previous     = t_current;
            v_previous     = v_current;
         }
      } while ( i<(N-1) ); 

      if(verbosity>=3) std::cout << "[NMRMath::CountZeroCrossings]: Done." << std::endl;
      // ClearVectors<T1,T2>(X,Y); 

      return cntr;   // return number of zero crossings  
   }

} // ::NMRMath 

#endif 
