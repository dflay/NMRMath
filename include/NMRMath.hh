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
      T GetMean(const int N,T *x){
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
      T GetRMS(const int N,T *x){
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
      T GetVariance(const int N,T *x){
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
      T GetStandardDeviation(const int N,T *x){
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
   template < typename T >
      T GetStandardErrorOfTheMean(const int N,T *x){
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
      T2 GetCovariance(const int N,T1 *x,T2 *y){
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
   template <typename T>
   int InitializeVector(size_t N,std::vector<T> &X){
      X.reserve(N);
      for(size_t i=0;i<N;i++) X.push_back(0);
      return 0;
   }
   //______________________________________________________________________________
   template <typename T>
   int ClearVector(std::vector<T> &X){
      const int N = X.size(); 
      for(int i=0;i<N;i++) X[i] = 0;
      return 0;
   }
   //______________________________________________________________________________
   int LeastSquaresFitting(std::vector<double> x,std::vector<double> y,double &intercept,double &slope,double &r){
      // linear regression to find slope b and y-intercept a of 
      // f(x) = intercept + slope*x 

      int rc=0;
      double num=0,rsq=0;

      const int N = x.size(); 
      double mean_x = GetMean<double>(x);
      double mean_y = GetMean<double>(y);
      double var_x  = GetVariance<double>(x);
      double var_y  = GetVariance<double>(y);
      double cov_xy = GetCovariance<double,double>(x,y);

      double ss_xx = ( (double)N )*var_x;
      double ss_yy = ( (double)N )*var_y;
      double ss_xy = ( (double)N )*cov_xy;

      double den = ss_xx*ss_yy;
      if(den==0){
	 // singular matrix. can't solve the problem.
	 intercept = 0;
	 slope     = 0;
	 r         = 0;
	 rc        = 1;
      }else{
	 slope     = cov_xy/var_x;
	 intercept = mean_y - slope*mean_x;
	 num       = ss_xy*ss_xy;
	 rsq       = num/den;
	 r         = sqrt(rsq);
      }

      return rc;
   }
   //______________________________________________________________________________
   int LeastSquaresFitting(const int N,std::vector<double> x,std::vector<double> y,double &intercept,double &slope,double &r){
      // we do this just in case we need to loop over N < x.size() 

      // linear regression to find slope b and y-intercept a of 
      // f(x) = a + bx 

      int rc=0;
      double num=0,rsq=0;

      double mean_x = GetMean<double>(N,x);
      double mean_y = GetMean<double>(N,y);
      double var_x  = GetVariance<double>(N,x);
      double var_y  = GetVariance<double>(N,y);
      double cov_xy = GetCovariance<double,double>(N,x,y);

      double ss_xx = ( (double)N )*var_x;
      double ss_yy = ( (double)N )*var_y;
      double ss_xy = ( (double)N )*cov_xy;

      double den = ss_xx*ss_yy;
      if(den==0){
	 // singular matrix. can't solve the problem.
	 intercept = 0;
	 slope     = 0;
	 r         = 0;
	 rc        = 1;
      }else{
	 slope     = cov_xy/var_x;
	 intercept = mean_y - slope*mean_x;
	 num       = ss_xy*ss_xy;
	 rsq       = num/den;
	 r         = sqrt(rsq);
      }

      return rc;
   }
   //______________________________________________________________________________
   int LeastSquaresFitting(const int N,double *x,double *y,double &intercept,double &slope,double &r){
      // we do this just in case we need to loop over N < x.size() 

      // linear regression to find slope b and y-intercept a of 
      // f(x) = a + bx 

      int rc=0;
      double num=0,rsq=0;

      double mean_x = GetMean<double>(N,x);
      double mean_y = GetMean<double>(N,y);
      double var_x  = GetVariance<double>(N,x);
      double var_y  = GetVariance<double>(N,y);
      double cov_xy = GetCovariance<double,double>(N,x,y);

      double ss_xx = ( (double)N )*var_x;
      double ss_yy = ( (double)N )*var_y;
      double ss_xy = ( (double)N )*cov_xy;

      double den = ss_xx*ss_yy;
      if(den==0){
	 // singular matrix. can't solve the problem.
	 intercept = 0;
	 slope     = 0;
	 r         = 0;
	 rc        = 1;
      }else{
	 slope     = cov_xy/var_x;
	 intercept = mean_y - slope*mean_x;
	 num       = ss_xy*ss_xy;
	 rsq       = num/den;
	 r         = sqrt(rsq);
      }

      return rc;
   }
   //______________________________________________________________________________
   double LinearInterpolationForTime(double y,double x0,double y0,double x1,double y1){
      // time can't be negative here, so we take the absolute values 
      // when calculating the slope  
      double b = std::abs(y-y0)/std::abs(y1-y0);
      double x = x0 + b*(x1-x0);
      return x;
   }
   //______________________________________________________________________________
   double LinearInterpolation(double x,double x0,double y0,double x1,double y1){
      double b = (x-x0)/(x1-x0);
      double y = y0 + b*(y1-y0); 
      return y;
   }
   //______________________________________________________________________________
   double GetTimeOfCrossing(int verbosity,int method,const int NN,
	 std::vector<double> X,std::vector<double> Y,
	 double t_current,double v_current,
	 double t_next   ,double v_next){
      // find the time of the zero crossing 
      // t_current,v_current, etc: (time,voltage) for midpoint method 
      // X,Y = (time,voltage) to use in the fitting procedure if doing linear interpolation or least squares 

      int ret_val=0;
      double t0=0;
      double v0=0,a=0,b=0,r=0;

      if (method==kMidpoint){
	 // method 1: take midpoint between t_current and t_next 
	 t0 = (t_current + t_next)/2.;
      } else if (method==kLinearInterpolation){
	 // method 2: get time at V = 0, linear interpolation  
	 t0 = LinearInterpolationForTime(v0,t_current,v_current,t_next,v_next);
      } else if(method==kLeastSquares){
	 // method 3: least squares fit to neighboring points
	 // to find fit parameters a and b in f(x) = a + bx 
	 ret_val = LeastSquaresFitting(NN,X,Y,a,b,r);
	 if (ret_val!=0) {
	    std::cout << "[NMRMath::GetTimeOfCrossing]: ERROR!  Least squares fit failed!" << std::endl;
	    t0 = -1;
	 } else {
	    if(b!=0){
	       t0 = -a/b;
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
	    t0 = LinearInterpolationForTime(v0,t_current,v_current,t_next,v_next);
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
   int StoreData(int verbosity,int i,int NPTS,
	 std::vector<double> time,std::vector<double> voltage,
	 std::vector<double> &X,std::vector<double> &Y){

      int N     = time.size();
      int start = i - NPTS/2;
      int end   = i + NPTS/2;

      double v_current = voltage[i];

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
      for(int j=start;j<=end;j++){
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
   int CountZeroCrossings(int verbosity,int method,double SampleFreq,double ExpFreq,
	 bool UseTimeRange,double tMin,double tMax,
	 std::vector<double> time,std::vector<double> voltage, 
	 std::vector<double> &tCross,std::vector<double> &vCross){
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

      double t0=0,t_current=0,t_previous=0,t_next=0;

      double v0=0,target=0;
      double v_prod=0,delta_v=0;
      double v_current=0,v_previous=0,v_next=0;

      int index=0;  // vector index
      std::vector<double> X,Y; 
      InitializeVector<double>(NPTS,X); 
      InitializeVector<double>(NPTS,Y); 

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
		  delta_v = std::abs(v_current-v_next); 
		  // fill vectors for fit method 
		  NPTSUseable = StoreData(verbosity,i,NPTS,time,voltage,X,Y); 
		  t0 = GetTimeOfCrossing(verbosity,method,NPTSUseable,X,Y,t_current,v_current,t_next,v_next); 
		  // fill vectors 
		  tCross[index] = t0; 
		  vCross[index] = v0; 
		  index++;  
	       }
	    } else {
	       // don't use the fit range
	       cntr++; // count the crossing 
	       delta_v = std::abs(v_current-v_next); 
	       // fill vectors for fit method 
	       NPTSUseable = StoreData(verbosity,i,NPTS,time,voltage,X,Y); 
	       t0 = GetTimeOfCrossing(verbosity,method,NPTSUseable,X,Y,t_current,v_current,t_next,v_next); 
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
	    // index          = 0; 
	    cntr_prev      = cntr; 
	    t_previous     = t_current;
	    v_previous     = v_current;
	 }
      } while ( i<(N-1) ); 

      if(verbosity>=3) std::cout << "[NMRMath::CountZeroCrossings]: Done." << std::endl;

      return cntr;   // return number of zero crossings  
   }
   //______________________________________________________________________________
   int ApplyOffset(double offset,std::vector<double> &voltage){
      // apply an offset to the data
      const int N = voltage.size(); 
      for(int i=0;i<N;i++) voltage[i] -= offset;
      return 0; 
   }
   //______________________________________________________________________________
   double GetTDiff(int nzc,std::vector<double> tCross,double &delta_t_even_nc,double &delta_t_odd_nc,bool &OffsetFail){
      delta_t_odd_nc=0;
      delta_t_even_nc=0;
      int counter_odd=0,counter_even=0;
      int NN = nzc;
      if(NN<2){
	 std::cout << "[NMRMath::GetTDiff]: NOT ENOUGH DATA POINTS! Number of data points: " << NN << std::endl;
	 OffsetFail = true;
	 return -1;
      }else{
	 for(int i=1;i<NN;i++){
	    if(i%2==0){
	       delta_t_odd_nc  += (tCross[i]-tCross[i-1]);   // we're going by even/odd Zc -- not Nc, so we interchange odd/even wrt Nc. 
	       counter_odd++;
	    }else if(i%2!=0){
	       delta_t_even_nc += (tCross[i]-tCross[i-1]);
	       counter_even++;
	    }
	    // std::cout << delta_t_odd_nc << "\t" << delta_t_even_nc << std::endl;
	 }
	 delta_t_odd_nc  /= ( (double)counter_odd );
	 delta_t_even_nc /= ( (double)counter_even);
      }
      double t_diff = delta_t_odd_nc-delta_t_even_nc;
      // std::cout << "dt_odd  = " << delta_t_odd_nc  << std::endl;
      // std::cout << "dt_even = " << delta_t_even_nc << std::endl;
      // std::cout << "dt_diff = " << t_diff << std::endl;
      return t_diff;
   }
   //______________________________________________________________________________
   int CheckOffset(double offset_old,double offset_new,double t_diff_old,double t_diff_new,double slope){

      int is_nan_t_diff_old = isnan(t_diff_old);
      int is_nan_t_diff_new = isnan(t_diff_new);
      int is_nan_offset_old = isnan(offset_old);
      int is_nan_offset_new = isnan(offset_new);
      int is_nan_slope      = isnan(slope);

      int rc = 0,rc_tot=0;
      if(is_nan_t_diff_old){
	 rc = 1;
	 rc_tot++;
      }
      if(is_nan_t_diff_new){
	 rc = 2;
	 rc_tot++;
      }
      if(is_nan_offset_old){
	 rc = 3;
	 rc_tot++;
      }
      if(is_nan_offset_new){
	 rc = 4;
	 rc_tot++;
      }
      if(is_nan_slope){
	 rc = 5;
	 rc_tot++;
      }

      // int dummy=0; 
      if(rc_tot>0){
	 printf("[NMRFileManager]: WARNING: One of the offset values below is NAN! \n");
	 printf("                  offset_old: %.7E \n",offset_old);
	 printf("                  offset_new: %.7E \n",offset_new);
	 printf("                  t_diff_old: %.7E \n",t_diff_old);
	 printf("                  t_diff_new: %.7E \n",t_diff_new);
	 printf("                  slope:      %.7E \n",slope     );
	 // printf("Enter any number to continue: ");
	 // cin  >> dummy; 
      }

      return rc_tot;
   }
   //______________________________________________________________________________
   int GetBaselineZC(int verbosity,double input_offset,double SampleFreq,double ExpFreq,std::vector<double> time,std::vector<double> voltage,double &output_offset){

      // find the baseline using zero crossings to determine if time between crossings is constant,
      // as is expected for pure sine waves

      // make a local copy for analysis 
      const int NN = time.size();
      std::vector<double> TIME,VOLTAGE; 
      for(int i=0;i<NN;i++){
         TIME.push_back(time[i]);
         VOLTAGE.push_back(voltage[i]);
      }
 
      double T_exp  = 1./ExpFreq;
      double N_exp  = T_exp*SampleFreq;       // number of points for one period 
      int step      = (int)( (1./16.)*N_exp );  // skip 1/8 of a period 
      int NPTS      = step/2;                  // use step/2 

      if(verbosity>=3) std::cout << "[NMRFileManager::GetOffsetZC]: Finding additional offset..." << std::endl;

      if(verbosity>=3){
         std::cout << "[NMRFileManager::GetOffsetZC]: Parameters: " << std::endl;
         std::cout << "                               Points in fit: " << NPTS << std::endl;
         std::cout << "                               step size:     " << step << std::endl;
      }

      // settings for counting zero crossings
      bool UseRange = false; 
      double tMin   = 0;    // in seconds 
      double tMax   = 1;    // in seconds                         
      int type      = kLeastSquares; 

      const int NC = 1E+5; 
      std::vector<double> tCross,vCross;
      InitializeVector<double>(NC,tCross);
      InitializeVector<double>(NC,vCross);

      std::cout << tCross.size() << std::endl;  
 
      int rc=0,counter=0; 

      double t_even=0,t_odd=0; 
      double err   = 1E-16; 

      double t_diff_abs=0,t_diff_abs_2=0; 

      // first calculation 
      int nzc = CountZeroCrossings(verbosity,type,SampleFreq,ExpFreq,UseRange,tMin,tMax,TIME,VOLTAGE,tCross,vCross);

      bool OffsetFail = false;
      double t_diff_old = GetTDiff(nzc,tCross,t_even,t_odd,OffsetFail); 

      if(OffsetFail) return 1;

      double offset_old  = input_offset; 
      double offset_new  = offset_old*(1. + 0.01);  // add 1%  

      if(offset_old==0) offset_new = 1E-6; 

      if(verbosity>=4){
         std::cout << "----------------------------------------------------------------" << std::endl;
         std::cout << "First calculation: " << std::endl;
         printf("offset_old = %.7E \n",offset_old);
         printf("t_diff_old = %.7E \n",t_diff_old);
         printf("offset_new = %.7E \n",offset_new);
         std::cout << "----------------------------------------------------------------" << std::endl;
      }
  
      ApplyOffset(offset_new,VOLTAGE);
      nzc = CountZeroCrossings(verbosity,type,SampleFreq,ExpFreq,UseRange,tMin,tMax,TIME,VOLTAGE,tCross,vCross);

      double t_diff_new = GetTDiff(nzc,tCross,t_even,t_odd,OffsetFail); 

      if(OffsetFail) return 1;
 
      // reset the analysis data
      TIME.clear();
      VOLTAGE.clear(); 
      for(int i=0;i<NN;i++){
         TIME.push_back(time[i]); 
         VOLTAGE.push_back(voltage[i]); 
      }

      double slope = (t_diff_new - t_diff_old)/(offset_new - offset_old);

      rc = CheckOffset(offset_old,offset_new,t_diff_old,t_diff_new,slope); 
      if(rc!=0) return rc;

      if(verbosity>=4){
         std::cout << "----------------------------------------------------------------" << std::endl;
         std::cout << "Second calculation: " << std::endl;
         printf("offset_old = %.7E \n",offset_old);
         printf("t_diff_old = %.7E \n",t_diff_old);
         printf("offset_new = %.7E \n",offset_new);
         printf("t_diff_new = %.7E \n",t_diff_new);
         printf("slope      = %.7E \n",slope)     ; 
         std::cout << "----------------------------------------------------------------" << std::endl;
      }

      offset_new = offset_old - t_diff_old/slope;
      double root_diff  = fabs(offset_new-offset_old);

      rc = CheckOffset(offset_old,offset_new,t_diff_old,t_diff_new,slope); 
      if(rc>0) return 1;

      if(t_diff_new==0) offset_new = offset_old;          // if new time difference is identically zero, we use old offset. 

      if(verbosity>4){
         std::cout << "trial     = 0"                   << std::endl; 
         printf("slope      = %.7E \n",slope); 
         printf("offset_old = %.7E   ",offset_old);
         printf("t_diff_old = %.7E \n",t_diff_old);
         printf("offset_new = %.7E \n",offset_new);
         printf("diff(offset-offset_prev) = %.7E \n",root_diff);
         std::cout << "----------------------------------------------------------------" << std::endl;
      }

      // we can't have a massive offset -- that's just not true. 
      if( fabs(offset_new)>0.005 ) offset_new = offset_old; 

      // update values 
      offset_old = offset_new; 
      t_diff_old = t_diff_new; 

      // clear arrays before starting
      // tCross.clear(); 
      // vCross.clear(); 

      int LIMIT = 20; 

      do{ 
         offset_new = offset_old - t_diff_old/slope;
         // check the new offset  
         rc = CheckOffset(offset_old,offset_new,t_diff_old,t_diff_new,slope); 
         if(rc>0) break; 
         ApplyOffset(offset_new,VOLTAGE); 
         nzc = CountZeroCrossings(verbosity,type,SampleFreq,ExpFreq,UseRange,tMin,tMax,TIME,VOLTAGE,tCross,vCross);
         t_diff_new = GetTDiff(nzc,tCross,t_even,t_odd,OffsetFail); 
         slope      = (t_diff_new - t_diff_old)/(offset_new - offset_old);
         root_diff  = fabs(offset_new - offset_old); 
         if(verbosity>4){ 
            std::cout << "trial     = " << counter+1 << std::endl; 
            printf("slope      = %.7E \n",slope); 
            printf("offset_old = %.7E   ",offset_old); 
            printf("t_diff_old = %.7E \n",t_diff_old);
            printf("offset_new = %.7E   ",offset_new); 
            printf("t_diff_new = %.7E \n",t_diff_new); 
            printf("diff(offset-offset_prev) = %.7E \n",root_diff);
            std::cout << "----------------------------------------------------------------" << std::endl;
         }
         if(OffsetFail) break;
         // set up for next calc 
         // reset the pulse (to try new offset on)   
         VOLTAGE.clear(); 
         VOLTAGE = voltage; 
         // for(int i=0;i<NN;i++) VOLTAGE.push_back(voltage[i]); 
         t_diff_abs   = fabs(t_diff_new); 
         t_diff_abs_2 = fabs(t_diff_new-t_diff_old); 
         t_diff_old   = t_diff_new;
         offset_old   = offset_new;
         counter++; 
         // clear vectors  
         // tCross.clear();
         // vCross.clear();
         std::cout << t_diff_abs << " " << err << std::endl;
      }while( (t_diff_abs>err)&&(counter<LIMIT) ); 

      if(counter==LIMIT){
         if(verbosity>=3){
            std::cout << "Reached max counter.  Giving up." << std::endl;
            printf("offset_old: %.7E \n",offset_old);
            printf("offset_new: %.7E \n",offset_new);
            printf("t_diff_old: %.7E \n",t_diff_old);
            printf("t_diff_new: %.7E \n",t_diff_new);
            printf("slope:      %.7E \n",slope     );
         }
      }

      if(verbosity>=3) std::cout << "[NMRFileManager::GetOffsetZC]: Done." << std::endl;

      TIME.clear();
      VOLTAGE.clear();
      // tCross.clear();
      // vCross.clear(); 

      if(rc>0 || OffsetFail){
         output_offset = 0;
      }else{
         output_offset = offset_new;
      }

      return rc; 
   }


} // ::NMRMath 

#endif 
