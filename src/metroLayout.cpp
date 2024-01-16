#include <Rcpp.h>
using namespace Rcpp;

double angle_between_edges(NumericVector pvec,NumericVector qvec) {
  if(pvec[0]==qvec[0] && pvec[1]==qvec[1]){
    return 0.0;
  }
  double dot_pq = pvec[0] * qvec[0] + pvec[1] * qvec[1];
  double mag_pq = sqrt(pvec[0]*pvec[0]+pvec[1]*pvec[1]) * sqrt(qvec[0]*qvec[0]+qvec[1]*qvec[1]);

  if(dot_pq/mag_pq < -0.99){
    return M_PI;
  }
  if(dot_pq/mag_pq > 0.99){
    return 0.0;
  }

  return acos(dot_pq/mag_pq);
}

// [[Rcpp::export]]
double criterion_angular_resolution(List adj,NumericMatrix xy){
  int n = adj.length();
  double crit=0;
  double elen;
  double angle;

  for(int i=0;i<n;++i){
    IntegerVector Ni = adj[i];
    if(Ni.length()!=1){
      NumericMatrix edges_xy(Ni.length(),4);
      for(int j=0;j<Ni.length();++j){
        edges_xy(j,0) = xy(i,0);
        edges_xy(j,1) = xy(i,1);
        edges_xy(j,2) = xy(Ni[j],0);
        edges_xy(j,3) = xy(Ni[j],1);

        edges_xy(j,2) -= edges_xy(j,0);
        edges_xy(j,3) -= edges_xy(j,1);

        elen = sqrt(edges_xy(j,2) * edges_xy(j,2) + edges_xy(j,3) * edges_xy(j,3));
        edges_xy(j,2) /= elen;
        edges_xy(j,3) /= elen;
      }
      for(int p=0;p<(Ni.length()-1);++p){
        for(int q=p+1;q<Ni.length();++q){
          NumericVector pvec = {edges_xy(p,2),edges_xy(p,3)};
          NumericVector qvec = {edges_xy(q,2),edges_xy(q,3)};
          angle = angle_between_edges(pvec,qvec);

          crit += std::abs(2.0 * M_PI / Ni.length() - angle);
        }
      }
    }
  }
  return crit;
}

// [[Rcpp::export]]
double criterion_edge_length(IntegerMatrix el, NumericMatrix xy, double lg){
  double crit = 0;
  double elen;
  NumericVector evec(4);
  for(int e=0;e<el.nrow();++e){
    evec = {xy(el(e,0),0),xy(el(e,0),1),
            xy(el(e,1),0),xy(el(e,1),1)};
    elen = sqrt((evec[0] - evec[2]) * (evec[0] - evec[2]) +
                (evec[1] - evec[3]) * (evec[1] - evec[3]));
    crit += std::abs(elen/lg - 1.0);
  }
  return crit;
}

// [[Rcpp::export]]
double criterion_balanced_edge_length(List adj_deg2,NumericMatrix xy){
  int n2 = adj_deg2.length();
  double crit = 0.0;

  if(n2==0){
    return 0.0;
  }
  NumericMatrix edges_xy(2,4);
  NumericMatrix elen(2);
  for(int i=0;i<n2;++i){
    IntegerVector Ni = adj_deg2[i];
    edges_xy(0,0) = xy(i,0);
    edges_xy(0,1) = xy(i,1);
    edges_xy(0,2) = xy(Ni[0],0);
    edges_xy(0,3) = xy(Ni[0],1);

    edges_xy(0,2) -= edges_xy(0,0);
    edges_xy(0,3) -= edges_xy(0,1);

    elen[0] = sqrt(edges_xy(0,2) * edges_xy(0,2) + edges_xy(0,3) * edges_xy(0,3));


    edges_xy(1,0) = xy(i,0);
    edges_xy(1,1) = xy(i,1);
    edges_xy(1,2) = xy(Ni[1],0);
    edges_xy(1,3) = xy(Ni[1],1);

    edges_xy(1,2) -= edges_xy(1,0);
    edges_xy(1,3) -= edges_xy(1,1);

    elen[1] = sqrt(edges_xy(1,2) * edges_xy(1,2) + edges_xy(1,3) * edges_xy(1,3));

    crit += std::abs(elen[1]-elen[0]);
  }
  return crit;
}

// [[Rcpp::export]]
double criterion_line_straightness(){
  return 0.0;
}

// [[Rcpp::export]]
double criterion_octilinearity(IntegerMatrix el,NumericMatrix xy){
  double crit = 0;
  NumericVector edge(4);
  for(int e=0;e<el.nrow();++e){
    edge = {xy(el(e,0),0),xy(el(e,0),1),xy(el(e,1),0),xy(el(e,1),1)};
    crit += std::abs(sin(4.0 * atan(std::abs(edge(1)-edge(3)) / std::abs(edge(0)-edge(2)))));
  }
  return crit;
}


double criterion_sum(List adj, IntegerMatrix el, List adj_deg2, NumericMatrix xy,
                     double lg, NumericVector w){
  double cn1 = criterion_angular_resolution(adj,xy);
  double cn2 = criterion_edge_length(el,xy,lg);
  double cn3 = criterion_balanced_edge_length(adj_deg2,xy);
  double cn4 = criterion_line_straightness();
  double cn5 = criterion_octilinearity(el,xy);

  return w[0] * cn1 + w[1] * cn2 + w[2] * cn3 + w[3] * cn4 + w[4] * cn5;
}

// [[Rcpp::export]]
NumericMatrix layout_as_metro_iter(List adj, IntegerMatrix el, List adj_deg2,
                                   NumericMatrix xy, NumericMatrix bbox,
                                   double l, double gr, NumericVector w, double bsize){
  double cur_crit, tmp_crit;
  double lg = l * gr;
  int n = adj.length();
  double x, y, xold, yold, xbest, ybest;
  bool running = true;
  int k = 0;
  NumericVector xoptions(8),yoptions(8);

  cur_crit = criterion_sum(adj,el,adj_deg2,xy,lg,w);

  while(running){
    running = false;
    for(int v=0;v<n;++v){
      xold = xy(v,0);
      yold = xy(v,1);
      xbest = xold;
      ybest = yold;
      xoptions = {xold-gr,xold,xold+gr,xold-gr,xold+gr,xold-gr,xold,xold+gr};
      yoptions = {yold+gr,yold+gr,yold+gr,yold,yold,yold-gr,yold-gr,yold-gr};
      for(int i=0;i<xoptions.length();++i){
        x = xoptions[i];
        y = yoptions[i];
        if((x>=bbox(v,0)) && (x<=bbox(v,2)) && (y>=bbox(v,1)) && (y<=bbox(v,3))){
          xy(v,0) = x;
          xy(v,1) = y;
          tmp_crit = criterion_sum(adj,el,adj_deg2,xy,lg,w);
          if(tmp_crit < cur_crit){
            cur_crit = tmp_crit;
            running = true;
            xbest = xy(v,0);
            ybest = xy(v,1);
          }
        }
      }
      xy(v,0) = xbest;
      xy(v,1) = ybest;
    }
    // Rcout << "run: " << k << " min: " << cur_crit << std::endl;
    k +=1;
  }
  return xy;
}
