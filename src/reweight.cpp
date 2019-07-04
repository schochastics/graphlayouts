#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector reweighting(IntegerMatrix el,List N_ranks) {
  int m = el.nrow();
  IntegerMatrix Nru;
  IntegerMatrix Nrv;
  IntegerVector Nru_vec;
  IntegerVector Nrv_vec;
  int u;
  int v;
  int ui;
  int vi;
  int maxu;
  int maxv;
  int maxi;
  NumericVector new_w(m);
  for(int e=0; e<m; ++e){
    u = el(e,0);
    v = el(e,1);
    Nru = as<IntegerMatrix>(N_ranks[u]);
    Nrv = as<IntegerMatrix>(N_ranks[v]);
    maxu = Rcpp::max(Nru(_, 1));
    maxv = Rcpp::max(Nrv(_, 1));
    maxi = max(IntegerVector::create(maxu,maxv));
    Nru_vec = Nru(_,0);
    Nrv_vec = Nrv(_,0);
    double jac_max = 0;
    for(int i=0;i<maxi;++i){
      if(i>maxu){
        ui = maxu;
      } else{
        ui = i;
      }
      if(i>maxv){
        vi = maxv;
      } else{
        vi = i;
      }
      IntegerVector Nru_sort = Nru_vec[seq_len(ui)];
      IntegerVector Nrv_sort = Nrv_vec[seq_len(vi)];
      Nru_sort = Nru_sort.sort();
      Nrv_sort = Nrv_sort.sort();
      int Nuv_int = intersect(Nru_sort,Nrv_sort).length();
      int Nuv_tot = union_(Nru_sort,Nrv_sort).length();
      double jac_temp = double(Nuv_int)/double(Nuv_tot);
      if(jac_temp>jac_max){
        jac_max = jac_temp;
      }
    }
    new_w[e]=jac_max;

  }
  return new_w;
}



