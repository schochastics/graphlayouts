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
    double stress_new=stress(xnew,W,D);
    double eps=(stress_old-stress_new)/stress_old;
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

