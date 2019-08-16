#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix sparseStress(NumericMatrix y,
                           NumericMatrix D,
                           List Rp,
                           IntegerVector pivots,
                           arma::sp_mat A,
                           int maxIter) {

  int j;

  double diff  = 1;
  double tdiff = 0;
  double xerr;
  double yerr;
  double denom;
  int iter=0;
  double tx;
  double ty;
  arma::rowvec deg;

  deg = sum(A,0);

  //number of nodes and pivots
  int n = y.nrow();
  int m = pivots.length();

  //initialize coordinates
  NumericMatrix x(n,2);
  for(int i=0;i<n;++i){
    x(i,0)=y(i,0);
    x(i,1)=y(i,1);
  }


  NumericMatrix W(n,m);
  NumericVector wsum(n);

  //reweighting
  for(int i=0;i<n;++i){
    for(int k=0;k<m;++k){
      int p = pivots[k];
      if(i!=p){
        double s=0;
        std::vector<int> Rpi = as<std::vector<int> >(Rp[k]);
        for(std::vector<int>::size_type l = 0; l!=Rpi.size(); ++l){
          j = Rpi[l];
          if(D(j,k)<=(D(i,k)/2)){
            s+=1;
          }
        }
        W(i,k) = s/(D(i,k)*D(i,k));

        if(A(i,p)!=1){
          wsum[i]+=W(i,k);
        }
      }
    }
  }

  while((diff > 0.0001) & (iter<maxIter)){
    iter+=1;
    diff=0;
    for(int i=0;i<n;i++){
      tx=0;
      ty=0;
      arma::sp_mat::const_col_iterator start = A.begin_col(i);
      arma::sp_mat::const_col_iterator end = A.end_col(i);
      for(arma::sp_mat::const_col_iterator k = start; k != end; ++k){
        j = k.row();
        denom = sqrt((x(i,0)-x(j,0))*(x(i,0)-x(j,0))+(x(i,1)-x(j,1))*(x(i,1)-x(j,1)));
        if(denom>0.001){
          tx+=x(j,0)+(x(i,0)-x(j,0))/denom;
          ty+=x(j,1)+(x(i,1)-x(j,1))/denom;
        }
      }

      for(int p=0;p<m;++p){
        j=pivots[p];
        if(A(i,j)==0){
          denom = sqrt((x(i,0)-x(j,0))*(x(i,0)-x(j,0))+(x(i,1)-x(j,1))*(x(i,1)-x(j,1)));
          if(denom>0.001){
            tx+=W(i,p)*(x(j,0)+(D(i,p)*(x(i,0)-x(j,0)))/denom);
            ty+=W(i,p)*(x(j,1)+(D(i,p)*(x(i,1)-x(j,1)))/denom);
          }
        }
      }

      tx=tx/(wsum[i]+deg[i]);
      ty=ty/(wsum[i]+deg[i]);

      if((tx > 0.1) & (x(i,0)>0.1)){
        xerr = std::abs((tx-x(i,0))/x(i,0));
      } else{
        xerr = 0;
      }
      if((ty>0.1) & (x(i,1)>0.1)){
        yerr = std::abs((ty-x(i,1))/x(i,1));
      } else{
        yerr = 0;
      }
      tdiff = xerr + yerr;

      diff+= tdiff;
      x(i,0) = tx;
      x(i,1) = ty;

    }

  }

  return x;
}

