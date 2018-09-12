#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double tolerance(NumericMatrix x,NumericMatrix y){
  int n = x.nrow();
  int dim = x.ncol();
  double eps=0;
  for(int i=0;i<n;++i){
      for(int d=0;d<dim;++d){
        eps+= sqrt((x(i,d)-y(i,d))*(x(i,d)-y(i,d)));
      }
  }
  eps=eps/(n*dim);
  return eps;
}

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
                     int dim,
                     int iter,
                     double tol) {
  int n = y.nrow();
  
  NumericMatrix x(n,dim);
  for(int i=0;i<n;++i){
    for(int d=0;d<dim;++d){
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
    NumericMatrix xnew(n,dim);
    for(int i=0;i<n;++i){
      for(int d=0; d<dim;++d){
        for(int j=0; j<n;++j){
          if(i!=j){
            double denom = sqrt((x(i,0)-x(j,0))*(x(i,0)-x(j,0))+(x(i,1)-x(j,1))*(x(i,1)-x(j,1)));
            xnew(i,d) += W(i,j)*(x(j,d)+D(i,j)*(x(i,d)-x(j,d))/denom);   
          }
        }
        xnew(i,d) = xnew(i,d)/wsum[i];
      }
    }
    double stress_new=stress(xnew,W,D);
    double eps=(stress_old-stress_new)/stress_old;
    if(eps<= tol){
      break;
    }
    stress_old=stress_new;
    for(int i=0;i<n;++i){
      for(int d=0;d<dim;++d){
        x(i,d)=xnew(i,d);
      }  
    }
  }
  return x;
}

