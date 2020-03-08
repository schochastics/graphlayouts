#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double stress(NumericMatrix x, NumericMatrix W, NumericMatrix D){
  double fct=0;
  int n=x.nrow();
  for(int i=0;i<(n-1);++i){
    for(int j=(i+1);j<n;++j){
      double denom = sqrt((x(i,0)-x(j,0))*(x(i,0)-x(j,0))+(x(i,1)-x(j,1))*(x(i,1)-x(j,1)));
      fct+=W(i,j)*(denom-D(i,j))*(denom-D(i,j));
    }
  }
  return fct;
}

// [[Rcpp::export]]
NumericMatrix stress_major(NumericMatrix y,
                     NumericMatrix W,
                     NumericMatrix D,
                     int iter,
                     double tol) {
  int n = y.nrow();

  NumericMatrix x(n,2);
  for(int i=0;i<n;++i){
    for(int d=0;d<2;++d){
      x(i,d)=y(i,d);
    }
  }

  NumericVector wsum(n);
  for(int i=0;i<n;++i){
    for(int j=0;j<n;++j){
      wsum[i]+=W(i,j);
    }
  }

  double stress_old = stress(x,W,D);

  for(int k=0; k<iter; ++k){
    NumericMatrix xnew(n,2); //out or in?
    for(int i=0;i<n;++i){
      for(int j=0; j<n;++j){
        if(i!=j){
          double denom = sqrt((x(i,0)-x(j,0))*(x(i,0)-x(j,0))+(x(i,1)-x(j,1))*(x(i,1)-x(j,1)));
          if(denom>0.00001){
            xnew(i,0) += W(i,j)*(x(j,0)+D(i,j)*(x(i,0)-x(j,0))/denom);
            xnew(i,1) += W(i,j)*(x(j,1)+D(i,j)*(x(i,1)-x(j,1))/denom);
          }
        }
      }
      xnew(i,0) = xnew(i,0)/wsum[i];
      xnew(i,1) = xnew(i,1)/wsum[i];
    }
    double stress_new = stress(xnew,W,D);
    double eps = (stress_old-stress_new)/stress_old;

    if(eps<= tol){
      break;
    }
    stress_old=stress_new;
    for(int i=0;i<n;++i){
      x(i,0)=xnew(i,0);
      x(i,1)=xnew(i,1);
    }
  }
  return x;
}

// [[Rcpp::export]]
NumericMatrix stress_radii(NumericMatrix y,
                           NumericMatrix W,
                           NumericMatrix D,
                           NumericVector r,
                           NumericVector tseq) {
  int n = y.nrow();
  int m = tseq.length();

  NumericMatrix x(n,2);
  for(int i=0;i<n;++i){
    for(int d=0;d<2;++d){
      x(i,d)=y(i,d);
    }
  }

  NumericVector wsum(n);
  for(int i=0;i<n;++i){
    for(int j=0;j<n;++j){
      wsum[i]+=W(i,j);
    }
  }

  NumericVector rpow(n);
  for(int i=0;i<n;++i){
    rpow[i] = 1/(r[i] * r[i]);
  }

  for(int s=0; s<m; ++s){
    double t = tseq[s];
    for(int k=0; k<1; ++k){
      NumericMatrix xnew(n,2); //out or in?
      for(int i=0;i<n;++i){
        for(int j=0; j<n;++j){
          if(i!=j){
            double bij = sqrt((x(i,0)-x(j,0))*(x(i,0)-x(j,0))+(x(i,1)-x(j,1))*(x(i,1)-x(j,1)));
            double ai = sqrt(x(i,0)*x(i,0)+x(i,1)*x(i,1));
            if(ai<0.0001){
              ai=0;
            } else{
              ai=1/ai;
            }
            if(bij<0.0001){
              bij=0;
            } else{
              bij=1/bij;
            }
            xnew(i,0) += (1-t) * W(i,j)*(x(j,0)+D(i,j)*(x(i,0)-x(j,0))*bij)+t*rpow[i]*(r[i]*x(i,0)*ai);
            xnew(i,1) += (1-t) * W(i,j)*(x(j,1)+D(i,j)*(x(i,1)-x(j,1))*bij)+t*rpow[i]*(r[i]*x(i,1)*ai);
          }
        }
        xnew(i,0) = xnew(i,0)/((1-t) * wsum[i] + t * rpow[i]);
        xnew(i,1) = xnew(i,1)/((1-t) * wsum[i] + t * rpow[i]);
      }
      for(int i=0;i<n;++i){
        x(i,0)=xnew(i,0);
        x(i,1)=xnew(i,1);
      }
    }
  }
  return x;
}

// [[Rcpp::export]]
NumericMatrix stress_focus(NumericMatrix y,
                           NumericMatrix W,
                           NumericMatrix D,
                           NumericMatrix Z,
                           NumericVector tseq,
                           int iter,
                           double tol) {
  int n = y.nrow();
  int m = tseq.length();

  NumericMatrix x(n,2);
  for(int i=0;i<n;++i){
    for(int d=0;d<2;++d){
      x(i,d)=y(i,d);
    }
  }

  NumericVector wsum(n);
  NumericVector zsum(n);
  for(int i=0;i<n;++i){
    for(int j=0;j<n;++j){
      wsum[i]+=W(i,j);
      zsum[i]+=Z(i,j);
    }
  }

  double stress_oldW = stress(x,W,D);
  // double stress_oldZ = stress(x,Z,D);
  double stress_old = stress_oldW;
  for(int s=0; s<m;++s){
    double t = tseq[s];

    for(int k=0; k<iter; ++k){
      NumericMatrix xnew(n,2); //out or in?
      for(int i=0;i<n;++i){
        for(int j=0; j<n;++j){
          if(i!=j){
            double denom = sqrt((x(i,0)-x(j,0))*(x(i,0)-x(j,0))+(x(i,1)-x(j,1))*(x(i,1)-x(j,1)));
            if(denom>0.00001){
              denom = 1/denom;
            } else{
              denom = 0;
            }
            xnew(i,0) += ((1-t)*W(i,j)+t*Z(i,j))*(x(j,0)+D(i,j)*(x(i,0)-x(j,0))*denom);
            xnew(i,1) += ((1-t)*W(i,j)+t*Z(i,j))*(x(j,1)+D(i,j)*(x(i,1)-x(j,1))*denom);
          }
        }
        xnew(i,0) = xnew(i,0)/((1-t)*wsum[i] + t*zsum[i]);
        xnew(i,1) = xnew(i,1)/((1-t)*wsum[i] + t*zsum[i]);
      }
      double stress_newW = stress(xnew,W,D);
      double stress_newZ = stress(xnew,Z,D);
      double stress_new = (1-t)*stress_newW + t*stress_newZ;

      for(int i=0;i<n;++i){
        x(i,0)=xnew(i,0);
        x(i,1)=xnew(i,1);
      }
      double eps=(stress_old-stress_new)/stress_old;
      if(eps<= tol){
        break;
      }
      stress_old=stress_new;
    }
  }
  return x;
}
