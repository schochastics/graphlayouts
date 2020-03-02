#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double constrained_stress(NumericMatrix x, NumericMatrix W, NumericMatrix D){
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
NumericMatrix constrained_stress_major(NumericMatrix y,
                           int dim,
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

  double stress_old = constrained_stress(x,W,D);

  for(int k=0; k<iter; ++k){
    NumericMatrix xnew(n,2); //out or in?
    for(int i=0;i<n;++i){
      if(dim==1){
        xnew(i,0) = y(i,0);
      } else{
        xnew(i,1) = y(i,1);
      }
      for(int j=0; j<n;++j){
        if(i!=j){
          double denom = sqrt((x(i,0)-x(j,0))*(x(i,0)-x(j,0))+(x(i,1)-x(j,1))*(x(i,1)-x(j,1)));
          if(denom>0.00001){
            if(dim==2){
              xnew(i,0) += W(i,j)*(x(j,0)+D(i,j)*(x(i,0)-x(j,0))/denom);
            } else{
              xnew(i,1) += W(i,j)*(x(j,1)+D(i,j)*(x(i,1)-x(j,1))/denom);
            }
          }
        }
      }
      if(dim==2){
        xnew(i,0) = xnew(i,0)/wsum[i];
      } else{
        xnew(i,1) = xnew(i,1)/wsum[i];
      }
    }
    double stress_new=constrained_stress(xnew,W,D);
    double eps=(stress_old-stress_new)/stress_old;
    if(eps<= tol){
      break;
    }
    stress_old=stress_new;
    for(int i=0;i<n;++i){
      x(i,0) = xnew(i,0);
      x(i,1) = xnew(i,1);
    }
  }
  return x;
}
