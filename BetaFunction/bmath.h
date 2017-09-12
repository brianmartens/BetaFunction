#ifndef BMATH_H_INCLUDED
#define BMATH_H_INCLUDED
#include "math.h"

using namespace std;

    const double ln_sqrt_2_pi = 0.91893853320467274178;
    const double g_pi = 3.14159265358979323846;
    const double num_e = 2.71828182845904523536;

  double abs(double x)
  {
      if (x >= 0)
        return x;
      else
        return -1*x;
  }
  double lanczos_ln_gamma(double z)
  {
    #define LG_g 7.0
    #define LG_N 9

  const double lct[LG_N+1] = {
     0.9999999999998099322768470047347,
    676.520368121885098567009190444019,
   -1259.13921672240287047156078755283,
    771.3234287776530788486528258894,
    -176.61502916214059906584551354,
     12.507343278686904814458936853,
    -0.13857109526572011689554707,
    9.984369578019570859563e-6,
    1.50563273514931155834e-7
  };
    double sum;
    double base;
    int i;
    if (z < 0.5) {
      // Use Euler's reflection formula:
      // Gamma(z) = Pi / [Sin[Pi*z] * Gamma[1-z]];
      return log(g_pi / sin(g_pi * z)) - lanczos_ln_gamma(1.0 - z);
    }
    z = z - 1.0;
    base = z + LG_g + 0.5;  // Base of the Lanczos exponential
    sum = 0;
    // We start with the terms that have the smallest coefficients and largest
    // denominator.
    for(i=LG_N; i>=1; i--) {
      sum += lct[i] / (z + ((double) i));
    }
    sum += lct[0];
    // Gamma[z] = Sqrt(2*Pi) * sum * base^[z + 0.5] / E^base
    return ((ln_sqrt_2_pi + log(sum)) - base) + log(base)*(z+0.5);
  }

  // Compute the Gamma function, which is e to the power of ln_gamma.
  double lanczos_gamma(double z)
  {
    return(exp(lanczos_ln_gamma(z)));
  }
  double low_incomplete_gamma(double s, double z)
  {
      double epsilon = .0000000001;
      int i = 1;
      double initial = (pow(z,s)*pow(num_e,-z))/s; // z^0 == 1
      for(;;)
      {
        double denom = s;
        for(int j=1 ; j<=i ; ++j)
        {
            denom = denom*(s+j);
        }
        double num = pow(z,s)*pow(num_e,-1*z)*pow(z,i);
        double test = num/denom;
        if(abs(test)<epsilon)
        {
            initial += test;
            break;
        }
        else
        {
            initial += test;
            ++i;
        }
      }
      return initial;
  }
  double regularized_gamma(double s, double z)
  {
      double inc = low_incomplete_gamma(s,z);
      double gam = lanczos_gamma(s);
      return inc/gam;
  }
  double error_function(double x)
  {
      return x < 0.0 ? -regularized_gamma(.5,x*x) : regularized_gamma(.5,x*x);
  }
  double Tdist(double t, int df)
  {
      double num = lanczos_gamma((df+1)/(double)2);
      double den = sqrt(g_pi*df)*lanczos_gamma(df/(double)2);
      double mult = (1 + pow(t,2)/(double)df);
      double exp =  (-1/(double)2)*(df+1);
      return (num/den)*(pow(mult,exp));
  }
  double Beta(double x, double y)
  {
      double gamx = lanczos_gamma(x);
      double gamy = lanczos_gamma(y);
      double gamxy = lanczos_gamma(x+y);
      return (gamx*gamy)/gamxy;
  }
  double i_beta_Lentz(double a, double b, double x)
  {
    #define MAXIT 100
    #define EPS 3.0e-7
    #define FPMIN 1.0e-30
    int m,m2;
    double aa,c,d,del,h,qab,qam,qap;
    qab=a+b;
    qap=a+1.0;
    qam=a-1.0;
    c=1.0;
    d=1.0-qab*x/qap;
    if (fabs(d) < FPMIN) d=FPMIN;
    d=1.0/d;
    h=d;
    for (m=1;m<=MAXIT;m++) {
        m2=2*m;
        aa=m*(b-m)*x/((qam+m2)*(a+m2));
        d=1.0+aa*d;
        if (fabs(d) < FPMIN) d=FPMIN;
        c=1.0+aa/c;
        if (fabs(c) < FPMIN) c=FPMIN;
        d=1.0/d;
        h *= d*c;
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
        d=1.0+aa*d;
        if (fabs(d) < FPMIN) d=FPMIN;
        c=1.0+aa/c;
        if (fabs(c) < FPMIN) c=FPMIN;
        d=1.0/d;
        del=d*c;
        h *= del;
        if (fabs(del-1.0) < EPS) break;
    }
    if (m > MAXIT)
    {
       cout << "a or b too big, or MAXIT too small in betacf";
       return 0.0;
    }
    else
    {
        return h;
    }
}
  double regularized_beta(double a,double b, double x)
  {
    double bt;
    if (x < 0.0 || x > 1.0) cout << "x out of range in regularized_beta" << endl;
    if (x == 0.0 || x == 1.0) bt=0.0;
    else
        bt=exp(lanczos_ln_gamma(a+b)-lanczos_ln_gamma(a)-lanczos_ln_gamma(b)+a*log(x)+b*log(1.0-x));
    if (x < (a+1.0)/(a+b+2.0)) // Use continued fraction directly.
        return bt*i_beta_Lentz(a,b,x)/a;
    else // Use continued fraction after making the symmetry transformation.
        return 1.0-bt*i_beta_Lentz(b,a,1.0-x)/b;
  }
  double t_dist_p(double t, double df)
  {
      double x = df/(pow(t,2)+df);
      double a = df/(double)2;
      double b = 1/(double)2;
      return double(1) - .5*regularized_beta(a,b,x);
  }
  double uw_pooled_z(double p1, int n1, double p2, int n2) // unweighted pooled z-test for two proportions (returns z value)
  {
    double z;
    double p_hat = ((p1*n1) + (p2*n2))/(n1+n2);
    double SE = sqrt(p_hat*(1.0-p_hat)*((1.0/n1) + (1.0/n2)));
    z = (p1-p2)/SE;
    return z;
  }
  double std_normal_CDF(double z)
  {
      return .5*(1.0 + error_function(z/sqrt(2.0)));
  }
  double f_dist_CDF(double v1, double v2, double f) // right tail probability
  {
      return regularized_beta(v2/(double)2, v1/(double)2, (v2/(v2+v1*f)));
  }
  double welch_uw_t_test(double xbar1, double xbar2, double var1, double var2, int n1, int n2) // unweighted welch's t-test
  {
      double meandiff = xbar1 - xbar2;
      double SE = sqrt((var1/n1) + (var2/n2));
      double t = meandiff/SE;
      double term1 = var1/n1;
      double term2 = var2/n2;
      double numerator = pow((term1 + term2),2);
      double term3 = pow(var1,2)/(pow(n1,2)*(n1-1));
      double term4 = pow(var2,2)/(pow(n2,2)*(n2-1));
      double denominator = term3 + term4;
      double df = numerator/denominator;
      return t_dist_p(t,df);
  }
  double chi_2_CDF(double df, double x)
  {
      return regularized_gamma(df/(double)2,x/(double)2);
  }
#endif // BMATH_H_INCLUDED
