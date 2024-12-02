#include <Rcpp.h>

using namespace Rcpp;

Function approx = Environment::namespace_env("stats")["approx"];

NumericVector appr(NumericVector t, NumericVector x, NumericVector new_t)
{
  List temp = approx(t, x, new_t);
  return(temp[1]);
}

int LCM_pair (int a, int b) {
  IntegerVector v = {a,b};
  int m = min(v), gcd = m, lcm;
  while (1) { if ((a % gcd == 0) & (b % gcd == 0)) {break;} else {gcd--;} }
  lcm = a*b/gcd;
  return lcm;
}

// [[Rcpp::export]]
NumericMatrix aCAPf1f2_cpp (List Jfunc1_ext_ls,
                            List Jfunc2_ext_ls,
                            double alpha,
                            int max_drv)
{
  NumericVector time = as<NumericMatrix>(Jfunc1_ext_ls[0])(_,0);
  int n_t = time.length();
  NumericMatrix h_func(n_t, 2);
  if (alpha == 0) {
    h_func(_,0) = h_func(_,1) = time;
  } else {
    NumericVector Jf1values, Jf2values;
    int i, j, ib, jb, i_dif, j_dif, optimal_i, optimal_j, k;
    double i_diff, j_diff, i_subb;
    List Jf1values_ls(max_drv), Jf2values_ls(max_drv);
    for (int b = 0; b<=max_drv-1; b++) {
      NumericVector temp1 = as<NumericMatrix>(Jfunc1_ext_ls[b])(_,1);
      NumericVector temp2 = as<NumericMatrix>(Jfunc2_ext_ls[b])(_,1);
      Jf1values_ls[b] = temp1;
      Jf2values_ls[b] = temp2;
    }

    IntegerMatrix lcm(max_drv + 1, max_drv + 1);
    for (i=1; i<=max_drv; i++) {
      for (j=1; j<=max_drv; j++) {
        if (i<=j) {lcm(i,j) = LCM_pair(i,j);} else {lcm(i,j) = lcm(j,i);}
      }
    }

    double cost_phase, cost_amp, sqrt_hprime, temp_cost, high_cost;
    int n_knots, i_sub, j_sub;
    NumericVector integrand;

    NumericMatrix partial_cost(n_t * n_t, (max_drv+1) * (max_drv+1));
    NumericMatrix partial_phase(max_drv+1, max_drv+1);
    for (i_dif = max_drv; i_dif >= 1; i_dif--) {
      i_diff = (double)(i_dif);
      for (j_dif = max_drv; j_dif >= 1; j_dif--) {
        j_diff = (double)(j_dif);
        sqrt_hprime = sqrt(j_diff/i_diff);
        cost_phase = sqrt_hprime * i_diff;
        partial_phase(i_dif, j_dif) = cost_phase;
        n_knots = lcm(i_dif,j_dif);
        i_sub = n_knots/i_dif; j_sub = n_knots/j_dif;
        i_subb = (double)(i_sub);
        Jf2values = Jf2values_ls[i_sub-1];
        Jf1values = Jf1values_ls[j_sub-1];
        for (j=1; j<=n_t-1; j++) {
          jb = j-j_dif;
          if (jb>=0) {
            for (i=1; i<=n_t-1; i++) {
              ib = i-i_dif;
              if (ib>=0) {
                integrand = sqrt_hprime *
                  Jf1values[seq(jb*j_sub, j*j_sub)] *
                  Jf2values[seq(ib*i_sub, i*i_sub)];
                cost_amp = sum(integrand[seq(0,n_knots-1)] + integrand[seq(1,n_knots)])/2 / i_subb;
                partial_cost(i*n_t+j, i_dif*(max_drv+1)+j_dif) = cost_amp*alpha + cost_phase*(1-alpha);
              }
            }
          }
        }
      }
    }

    NumericMatrix mm(2,1); NumericVector vv(2);
    List path(n_t*n_t);
    for (j=0;j<=n_t-1;j++) {path[n_t*j] = path[j] = mm;}

    NumericMatrix cost(n_t, n_t);
    cost(_,0) = cost(0,_) = rep(-99999,n_t); cost(0,0) = 0;
    if (alpha != 1) {
      for (i=1; i<=n_t-1; i++) {
        for (j=1; j<=n_t-2; j++) {
          if (i==(n_t-1)) {j=n_t-1;}
          k = 0;
          for (i_dif=max_drv; i_dif>=1; i_dif--) {
            ib = i-i_dif;
            if (ib>=0) {
              for (j_dif=max_drv; j_dif>=1; j_dif--) {
                jb = j-j_dif;
                if (jb>=0) {
                  temp_cost = cost(ib,jb) + partial_cost(i*n_t+j, i_dif*(max_drv+1)+j_dif);
                  if (k==0) {
                    optimal_i = ib; optimal_j = jb; high_cost = temp_cost; k = 1;
                  } else {
                    if (temp_cost > high_cost) {
                      optimal_i = ib; optimal_j = jb; high_cost = temp_cost;
                    }
                  }
                }
              }
            }
          }
          vv = {(double)(i),(double)(j)}; cost(i,j) = high_cost;
          path[i*n_t+j] = cbind(as<NumericMatrix>(path[optimal_i*n_t+optimal_j]),vv);
        }
      }
    } else {
      NumericMatrix cost_Ponly(n_t, n_t);
      cost_Ponly(_,0) = cost_Ponly(0,_) = rep(-99999,n_t); cost_Ponly(0,0) = 0;
      double temp_cost_Ponly, high_cost_Ponly;
      for (i=1; i<=n_t-1; i++) {
        for (j=1; j<=n_t-2; j++) {
          if (i==(n_t-1)) {j=n_t-1;}
          k=0;
          for (i_dif=max_drv; i_dif>=1; i_dif--) {
            ib = i-i_dif;
            if (ib>=0) {
              for (j_dif=max_drv; j_dif>=1; j_dif--) {
                jb = j-j_dif;
                if (jb>=0) {
                  temp_cost = cost(ib,jb) + partial_cost(i*n_t+j, i_dif*(max_drv+1)+j_dif);
                  temp_cost_Ponly = cost_Ponly(ib,jb) + partial_phase(i_dif, j_dif);
                  if (k==0) {
                    optimal_i = ib; optimal_j = jb; high_cost = temp_cost; high_cost_Ponly = temp_cost_Ponly;
                    k = 1;
                  } else {
                    if (temp_cost > high_cost) {
                      optimal_i = ib; optimal_j = jb; high_cost = temp_cost; high_cost_Ponly = temp_cost_Ponly;
                    } else if (temp_cost == high_cost && temp_cost_Ponly > high_cost_Ponly) {
                      optimal_i = ib; optimal_j = jb; high_cost = temp_cost; high_cost_Ponly = temp_cost_Ponly;
                    }
                  }
                }
              }
            }
          }
          vv = {(double)(i),(double)(j)}; cost(i,j) = high_cost; cost_Ponly(i,j) = high_cost_Ponly;
          path[i*n_t+j] = cbind(as<NumericMatrix>(path[optimal_i*n_t+optimal_j]),vv);
        }
      }
    }

    NumericMatrix final_path = path[n_t*n_t-1];
    final_path = final_path/(n_t-1);
    h_func(_,0) = time; h_func(_,1) = as<NumericVector>(appr(final_path(0,_),final_path(1,_),time));
  }

  return (h_func);
}
