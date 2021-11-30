#include <RcppArmadillo.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace arma;
using namespace Eigen;

SEXP pmvnorm(NumericVector lb, 
             NumericVector ub, 
             NumericVector m, 
             const arma::mat& s){
  Environment pkg = Environment::namespace_env("mvtnorm");
  Function f = pkg["pmvnorm"];
  return f(Named("lower") = lb, 
           _["upper"] = ub, 
           _["mean"] = m, 
           _["sigma"] = s);
}
// [[Rcpp::export(.CalAPCEj_ordinal_rcpp)]]
List CalAPCEj_ordinal_rcpp(const arma::mat& X,
                          const arma::mat& ZX,
                          double rho,
                          const arma::vec& beta,
                          const arma::vec& alpha,
                          NumericVector theta,
                          NumericVector delta,
                          int p,
                          int k,
                          int n,
                          NumericVector s1,
                          NumericVector s2,
                          NumericVector s3,
                          NumericVector s4,
                          NumericVector s5,
                          int lZX,
                          const arma::mat& C,
                          const arma::mat& C_binary,
                          int only_optimal_decision,
                          const arma::urowvec& D,
                          const arma::urowvec& D_binary,
                          arma::urowvec& dmf,
                          int dmf_null,
                          int only_dmf_fair){
  NumericVector theta1_extend = theta[Range(0,k-1)];
  NumericVector theta0_extend = theta[Range(k,2*k-1)];
  theta1_extend.push_front(R_NegInf);
  theta1_extend.push_back(R_PosInf);
  theta0_extend.push_front(R_NegInf);
  theta0_extend.push_back(R_PosInf);
  
  // E[e_r(X_i)]: calculate  Pr(R=r) for r=0,...,k+1 for the overall and subgroups defined by gender (1,0) and race (white=1)
  mat pj_R (5, k+2);
  mat g_d (n, k+1);
  g_d.fill(0);
  mat g_d_binary (n, 2);
  g_d_binary.fill(0);
  vec mu_R = X*alpha;
  arma::cube dmfj_DR;
  dmfj_DR.set_size(5,k+1,k+2);
  for (int r = 0; r<k; r++) {
    rowvec w = C.row(r+1);
    NumericVector v (n);
    for (int i = 0; i<n; i++) {
      v[i] = R::pnorm(delta[r+1], mu_R[i], 1.0, 1, 0) - R::pnorm(delta[r], mu_R[i], 1.0, 1, 0);
    }
    g_d += as<colvec>(v) * w;
    rowvec w_binary = C_binary.row(r+1);
    g_d_binary += as<colvec>(v) * w_binary;
    NumericVector v1 = v[s1];
    NumericVector v2 = v[s2];
    NumericVector v3 = v[s3];
    NumericVector v4 = v[s4];
    NumericVector v5 = v[s5];
    pj_R(0,r+1) = mean(v1);
    pj_R(1,r+1) = mean(v2);
    pj_R(2,r+1) = mean(v3);
    pj_R(3,r+1) = mean(v4);
    pj_R(4,r+1) = mean(v5);
    for (int d = 0; d<(k+1); d++) {
      NumericVector d_tilde (n);
      for (int i = 0; i<n; i++) {
        if(dmf[i] == d) {
          d_tilde[i] = 1;
        } else {
          d_tilde[i] = 0;
        }
      }
      NumericVector d_tilde1 = d_tilde[s1];
      NumericVector d_tilde2 = d_tilde[s2];
      NumericVector d_tilde3 = d_tilde[s3];
      NumericVector d_tilde4 = d_tilde[s4];
      NumericVector d_tilde5 = d_tilde[s5];
      NumericVector dmfj_DR1 = v1 * d_tilde1;
      NumericVector dmfj_DR2 = v2 * d_tilde2;
      NumericVector dmfj_DR3 = v3 * d_tilde3;
      NumericVector dmfj_DR4 = v4 * d_tilde4;
      NumericVector dmfj_DR5 = v5 * d_tilde5;
      dmfj_DR(0,d,r+1) = sum(dmfj_DR1)/sum(v1);
      dmfj_DR(1,d,r+1) = sum(dmfj_DR2)/sum(v2);
      dmfj_DR(2,d,r+1) = sum(dmfj_DR3)/sum(v3);
      dmfj_DR(3,d,r+1) = sum(dmfj_DR4)/sum(v4);
      dmfj_DR(4,d,r+1) = sum(dmfj_DR5)/sum(v5);
    }
  }
  NumericVector v0 (n);
  rowvec w0 = C.row(0);
  for (int i = 0; i<n; i++) {
    v0[i] = R::pnorm(delta[0], mu_R[i], 1.0, 1, 0);
  }
  g_d += as<colvec>(v0) * w0;
  rowvec w0_binary = C_binary.row(0);
  g_d_binary += as<colvec>(v0) * w0_binary;
  NumericVector v01 = v0[s1];
  NumericVector v02 = v0[s2];
  NumericVector v03 = v0[s3];
  NumericVector v04 = v0[s4];
  NumericVector v05 = v0[s5];
  pj_R(0,0) = mean(v01);
  pj_R(1,0) = mean(v02);
  pj_R(2,0) = mean(v03);
  pj_R(3,0) = mean(v04);
  pj_R(4,0) = mean(v05);
  
  NumericVector vL (n);
  rowvec wL = C.row(k+1);
  for (int i = 0; i<n; i++) {
    vL[i] = (1 - R::pnorm(delta[k], mu_R[i], 1.0, 1, 0));
  }
  g_d += as<colvec>(vL) * wL;
  rowvec wL_binary = C_binary.row(k+1);
  g_d_binary += as<colvec>(vL) * wL_binary;
  NumericVector vL1 = vL[s1];
  NumericVector vL2 = vL[s2];
  NumericVector vL3 = vL[s3];
  NumericVector vL4 = vL[s4];
  NumericVector vL5 = vL[s5];
  pj_R(0,k+1) = mean(vL1);
  pj_R(1,k+1) = mean(vL2);
  pj_R(2,k+1) = mean(vL3);
  pj_R(3,k+1) = mean(vL4);
  pj_R(4,k+1) = mean(vL5);
  
  for (int d = 0; d<(k+1); d++) {
    NumericVector d_tilde0L (n);
    for (int i = 0; i<n; i++) {
      if(dmf[i] == d) {
        d_tilde0L[i] = 1;
      } else {
        d_tilde0L[i] = 0;
      }
    }
    NumericVector d_tilde0L1 = d_tilde0L[s1];
    NumericVector d_tilde0L2 = d_tilde0L[s2];
    NumericVector d_tilde0L3 = d_tilde0L[s3];
    NumericVector d_tilde0L4 = d_tilde0L[s4];
    NumericVector d_tilde0L5 = d_tilde0L[s5];
    NumericVector dmfj_DR01 = v01 * d_tilde0L1;
    NumericVector dmfj_DR02 = v02 * d_tilde0L2;
    NumericVector dmfj_DR03 = v03 * d_tilde0L3;
    NumericVector dmfj_DR04 = v04 * d_tilde0L4;
    NumericVector dmfj_DR05 = v05 * d_tilde0L5;
    dmfj_DR(0,d,0) = sum(dmfj_DR01)/sum(v01);
    dmfj_DR(1,d,0) = sum(dmfj_DR02)/sum(v02);
    dmfj_DR(2,d,0) = sum(dmfj_DR03)/sum(v03);
    dmfj_DR(3,d,0) = sum(dmfj_DR04)/sum(v04);
    dmfj_DR(4,d,0) = sum(dmfj_DR05)/sum(v05);
    NumericVector dmfj_DRL1 = vL1 * d_tilde0L1;
    NumericVector dmfj_DRL2 = vL2 * d_tilde0L2;
    NumericVector dmfj_DRL3 = vL3 * d_tilde0L3;
    NumericVector dmfj_DRL4 = vL4 * d_tilde0L4;
    NumericVector dmfj_DRL5 = vL5 * d_tilde0L5;
    dmfj_DR(0,d,k+1) = sum(dmfj_DRL1)/sum(vL1);
    dmfj_DR(1,d,k+1) = sum(dmfj_DRL2)/sum(vL2);
    dmfj_DR(2,d,k+1) = sum(dmfj_DRL3)/sum(vL3);
    dmfj_DR(3,d,k+1) = sum(dmfj_DRL4)/sum(vL4);
    dmfj_DR(4,d,k+1) = sum(dmfj_DRL5)/sum(vL5);
  }
  // Optimal decision: row-wise maximum of g_d matrix
  urowvec g = index_max(g_d.t(), 0);
  mat delta_star_i(n, k+1);
  delta_star_i.fill(0);
  for (int i = 0; i<n; i++) {
    delta_star_i(i, g[i]) += 1;
  }
  vec g_d_i(n);
  g_d_i.fill(0);
  for (int i = 0; i<n; i++) {
    g_d_i(i) += g_d_binary(i, D_binary[i]);
  }
  if (dmf_null == 1) {
    for (int i = 0; i<n; i++) {
      dmf(i) += D_binary(i);
    }
  }
  vec g_dmf_i(n);
  g_dmf_i.fill(0);
  for (int i = 0; i<n; i++) {
    g_dmf_i(i) += g_d_binary(i, dmf[i]);
  }
  mat delta_star(5, k+1);
  delta_star.fill(0);
  int len1 = s1.length();
  int len2 = s2.length();
  int len3 = s3.length();
  int len4 = s4.length();
  int len5 = s5.length();
  for (int i = 0; i<len1; i++) {
    delta_star(0, g[s1[i]]) += 1/double(len1);
  }
  for (int i = 0; i<len2; i++) {
    delta_star(1, g[s2[i]]) += 1/double(len2);
  }
  for (int i = 0; i<len3; i++) {
    delta_star(2, g[s3[i]]) += 1/double(len3);
  }
  for (int i = 0; i<len4; i++) {
    delta_star(3, g[s4[i]]) += 1/double(len4);
  }
  for (int i = 0; i<len5; i++) {
    delta_star(4, g[s5[i]]) += 1/double(len5);
  }
  List out;
  if (only_optimal_decision==1) {
    // save results
    out = List::create(Named("Pj_R") = pj_R,
                            _["Optimal_D"] = delta_star,
                            _["Optimal_D_ind"] = delta_star_i,
                            _["Utility_g_d_ind"] = g_d_i,
                            _["Utility_g_dmf_ind"] = g_dmf_i);
  } else if (only_dmf_fair==1) {
    // save results
    out = List::create(Named("Pj_R") = pj_R,
                       _["Pj_dmf"] = dmfj_DR);
  } else {
    // E[Pr(D_i(z)=d|R_i=r,X_i)*e_r(X_i)]/E[e_r(X_i)] for all units (z=0,1 and d=0,1,...,k and r=0,1,...,k)
    arma::cube p1j_DR;
    arma::cube p0j_DR;
    p1j_DR.set_size(5,k+1,k+2);
    p0j_DR.set_size(5,k+1,k+2);
    NumericVector E1 (n);
    NumericVector E0 (n);
    E1.fill(0);
    E0.fill(0);
    vec b0 (p);
    vec b1 (lZX);
    vec b2 (p);
    b0.zeros();
    b1.zeros();
    b2.zeros();
    for (int i = 1; i<(lZX+1); i++) {
      b1[i-1] = beta[p+i];
    }
    for (int i = 1; i<(p+1); i++) {
      b0[i-1] = beta[i];
      b2[i-1] = beta[i];
    }
    if (rho==0){
      vec mu1_D = beta[0]+ X*b0 + ZX*b1;
      vec mu1_R = X*alpha;
      vec mu2_D = X*b2;
      vec mu2_R = mu1_R;
      for (int d = -1; d<k; d++) {
        for (int r = 0; r<k; r++) {
          // Pr(D_i(1)=d|R_i=r,X_i)*Pr(R_i=r|X_i) = Pr(D_i(1)=d|R_i=r,X_i)*e_r(X_i)
          // Pr(D_i(0)=d|R_i=r,X_i)*Pr(R_i=r|X_i) = Pr(D_i(0)=d|R_i=r,X_i)*e_r(X_i)
          NumericVector ind_1j_DR (n);
          NumericVector ind_0j_DR (n);
          for (int i = 0; i<n; i++) {
            ind_1j_DR[i] = (R::pnorm(theta1_extend[d+2], mu1_D[i], 1.0, 1, 0)-R::pnorm(theta1_extend[d+1], mu1_D[i], 1.0, 1, 0))*(R::pnorm(delta[r+1], mu1_R[i], 1.0, 1, 0)-R::pnorm(delta[r], mu1_R[i], 1.0, 1, 0));
            ind_0j_DR[i] = (R::pnorm(theta0_extend[d+2], mu2_D[i], 1.0, 1, 0)-R::pnorm(theta0_extend[d+1], mu2_D[i], 1.0, 1, 0))*(R::pnorm(delta[r+1], mu2_R[i], 1.0, 1, 0)-R::pnorm(delta[r], mu2_R[i], 1.0, 1, 0));
          }
          NumericVector ind_1j_DR1 = ind_1j_DR[s1];
          NumericVector ind_1j_DR2 = ind_1j_DR[s2];
          NumericVector ind_1j_DR3 = ind_1j_DR[s3];
          NumericVector ind_1j_DR4 = ind_1j_DR[s4];
          NumericVector ind_1j_DR5 = ind_1j_DR[s5];
          p1j_DR(0,d+1,r+1) = mean(ind_1j_DR1);
          p1j_DR(1,d+1,r+1) = mean(ind_1j_DR2);
          p1j_DR(2,d+1,r+1) = mean(ind_1j_DR3);
          p1j_DR(3,d+1,r+1) = mean(ind_1j_DR4);
          p1j_DR(4,d+1,r+1) = mean(ind_1j_DR5);
          NumericVector ind_0j_DR1 = ind_0j_DR[s1];
          NumericVector ind_0j_DR2 = ind_0j_DR[s2];
          NumericVector ind_0j_DR3 = ind_0j_DR[s3];
          NumericVector ind_0j_DR4 = ind_0j_DR[s4];
          NumericVector ind_0j_DR5 = ind_0j_DR[s5];
          p0j_DR(0,d+1,r+1) = mean(ind_0j_DR1);
          p0j_DR(1,d+1,r+1) = mean(ind_0j_DR2);
          p0j_DR(2,d+1,r+1) = mean(ind_0j_DR3);
          p0j_DR(3,d+1,r+1) = mean(ind_0j_DR4);
          p0j_DR(4,d+1,r+1) = mean(ind_0j_DR5);
          if (d==r) {
            E1 = E1 + ind_1j_DR;
            E0 = E0 + ind_0j_DR;
          }
        }
        NumericVector ind_1j_DR0 (n);
        NumericVector ind_0j_DR0 (n);
        for (int i = 0; i<n; i++) {
          ind_1j_DR0[i] = (R::pnorm(theta1_extend[d+2], mu1_D[i], 1.0, 1, 0)-R::pnorm(theta1_extend[d+1], mu1_D[i], 1.0, 1, 0))*(R::pnorm(delta[0], mu1_R[i], 1.0, 1, 0));
          ind_0j_DR0[i] = (R::pnorm(theta0_extend[d+2], mu2_D[i], 1.0, 1, 0)-R::pnorm(theta0_extend[d+1], mu2_D[i], 1.0, 1, 0))*(R::pnorm(delta[0], mu2_R[i], 1.0, 1, 0));
        }
        NumericVector ind_1j_DR01 = ind_1j_DR0[s1];
        NumericVector ind_1j_DR02 = ind_1j_DR0[s2];
        NumericVector ind_1j_DR03 = ind_1j_DR0[s3];
        NumericVector ind_1j_DR04 = ind_1j_DR0[s4];
        NumericVector ind_1j_DR05 = ind_1j_DR0[s5];
        p1j_DR(0,d+1,0) = mean(ind_1j_DR01);
        p1j_DR(1,d+1,0) = mean(ind_1j_DR02);
        p1j_DR(2,d+1,0) = mean(ind_1j_DR03);
        p1j_DR(3,d+1,0) = mean(ind_1j_DR04);
        p1j_DR(4,d+1,0) = mean(ind_1j_DR05);
        NumericVector ind_0j_DR01 = ind_0j_DR0[s1];
        NumericVector ind_0j_DR02 = ind_0j_DR0[s2];
        NumericVector ind_0j_DR03 = ind_0j_DR0[s3];
        NumericVector ind_0j_DR04 = ind_0j_DR0[s4];
        NumericVector ind_0j_DR05 = ind_0j_DR0[s5];
        p0j_DR(0,d+1,0) = mean(ind_0j_DR01);
        p0j_DR(1,d+1,0) = mean(ind_0j_DR02);
        p0j_DR(2,d+1,0) = mean(ind_0j_DR03);
        p0j_DR(3,d+1,0) = mean(ind_0j_DR04);
        p0j_DR(4,d+1,0) = mean(ind_0j_DR05);
        if ((d+1)==0) {
          E1 = E1 + ind_1j_DR0;
          E0 = E0 + ind_0j_DR0;
        }
        NumericVector ind_1j_DRL (n);
        NumericVector ind_0j_DRL (n);
        for (int i = 0; i<n; i++) {
          ind_1j_DRL[i] = (R::pnorm(theta1_extend[d+2], mu1_D[i], 1.0, 1, 0)-R::pnorm(theta1_extend[d+1], mu1_D[i], 1.0, 1, 0))*(1-R::pnorm(delta[k], mu1_R[i], 1.0, 1, 0));
          ind_0j_DRL[i] = (R::pnorm(theta0_extend[d+2], mu2_D[i], 1.0, 1, 0)-R::pnorm(theta0_extend[d+1], mu2_D[i], 1.0, 1, 0))*(1-R::pnorm(delta[k], mu2_R[i], 1.0, 1, 0));
        }
        NumericVector ind_1j_DRL1 = ind_1j_DRL[s1];
        NumericVector ind_1j_DRL2 = ind_1j_DRL[s2];
        NumericVector ind_1j_DRL3 = ind_1j_DRL[s3];
        NumericVector ind_1j_DRL4 = ind_1j_DRL[s4];
        NumericVector ind_1j_DRL5 = ind_1j_DRL[s5];
        p1j_DR(0,d+1,k+1) = mean(ind_1j_DRL1);
        p1j_DR(1,d+1,k+1) = mean(ind_1j_DRL2);
        p1j_DR(2,d+1,k+1) = mean(ind_1j_DRL3);
        p1j_DR(3,d+1,k+1) = mean(ind_1j_DRL4);
        p1j_DR(4,d+1,k+1) = mean(ind_1j_DRL5);
        NumericVector ind_0j_DRL1 = ind_0j_DRL[s1];
        NumericVector ind_0j_DRL2 = ind_0j_DRL[s2];
        NumericVector ind_0j_DRL3 = ind_0j_DRL[s3];
        NumericVector ind_0j_DRL4 = ind_0j_DRL[s4];
        NumericVector ind_0j_DRL5 = ind_0j_DRL[s5];
        p0j_DR(0,d+1,k+1) = mean(ind_0j_DRL1);
        p0j_DR(1,d+1,k+1) = mean(ind_0j_DRL2);
        p0j_DR(2,d+1,k+1) = mean(ind_0j_DRL3);
        p0j_DR(3,d+1,k+1) = mean(ind_0j_DRL4);
        p0j_DR(4,d+1,k+1) = mean(ind_0j_DRL5);
      }
    } else {
      mat cov_mat (2,2);
      cov_mat.fill(rho);
      cov_mat.diag() += (1-rho);
      arma::cube p1ij_DR;
      arma::cube p0ij_DR;
      p1ij_DR.set_size(k+1,k+2,n);
      p0ij_DR.set_size(k+1,k+2,n);
      for (int i = 0; i<n; i++) {
        rowvec x = X.row(i);
        rowvec zx = ZX.row(i);
        double mu1i_D = beta[0] + dot(x, b0) + dot(zx, b1);
        double mu1i_R = dot(x, alpha);
        double mu2i_D = dot(x, b2);
        double mu2i_R = mu1i_R;
        for (int d = -1; d<k; d++) {
          for (int r = 0; r<k; r++) {
            p1ij_DR(d+1,r+1,i) = as<double>(pmvnorm(NumericVector::create(theta1_extend[d+1],delta[r]), 
                                            NumericVector::create(theta1_extend[d+2], delta[r+1]), 
                                            NumericVector::create(mu1i_D,mu1i_R), 
                                            cov_mat));
            p0ij_DR(d+1,r+1,i) = as<double>(pmvnorm(NumericVector::create(theta0_extend[d+1],delta[r]), 
                                            NumericVector::create(theta0_extend[d+2], delta[r+1]), 
                                            NumericVector::create(mu2i_D,mu2i_R), 
                                            cov_mat));
          }
          p1ij_DR(d+1,0,i) = as<double>(pmvnorm(NumericVector::create(theta1_extend[d+1],R_NegInf), 
                                        NumericVector::create(theta1_extend[d+2], delta[0]), 
                                        NumericVector::create(mu1i_D,mu1i_R), 
                                        cov_mat));
          p0ij_DR(d+1,0,i) = as<double>(pmvnorm(NumericVector::create(theta0_extend[d+1],R_NegInf), 
                                        NumericVector::create(theta0_extend[d+2], delta[0]), 
                                        NumericVector::create(mu2i_D,mu2i_R), 
                                        cov_mat));
          p1ij_DR(d+1,k+1,i) = as<double>(pmvnorm(NumericVector::create(theta1_extend[d+1],delta[k]), 
                                        NumericVector::create(theta1_extend[d+2], R_PosInf), 
                                        NumericVector::create(mu1i_D,mu1i_R), 
                                        cov_mat));
          p0ij_DR(d+1,k+1,i) = as<double>(pmvnorm(NumericVector::create(theta0_extend[d+1],delta[k]), 
                                        NumericVector::create(theta0_extend[d+2], R_PosInf), 
                                        NumericVector::create(mu2i_D,mu2i_R), 
                                        cov_mat));
        }
        mat e1 = p1ij_DR.slice(i);
        mat e0 = p0ij_DR.slice(i);
        E1[i] = sum(e1.diag());
        E0[i] = sum(e0.diag());
      }
      arma::cube c = p1ij_DR.slices(as<uvec>(s1));
      arma::cube cube_mean = mean(c, 2);
      p1j_DR.row(0) = cube_mean.slice(0);
      c = p1ij_DR.slices(as<uvec>(s2));
      cube_mean = mean(c, 2);
      p1j_DR.row(1) = cube_mean.slice(0);
      c = p1ij_DR.slices(as<uvec>(s3));
      cube_mean = mean(c, 2);
      p1j_DR.row(2) = cube_mean.slice(0);
      c = p1ij_DR.slices(as<uvec>(s4));
      cube_mean = mean(c, 2);
      p1j_DR.row(3) = cube_mean.slice(0);
      c = p1ij_DR.slices(as<uvec>(s5));
      cube_mean = mean(c, 2);
      p1j_DR.row(4) = cube_mean.slice(0);
      
      c = p0ij_DR.slices(as<uvec>(s1));
      cube_mean = mean(c, 2);
      p0j_DR.row(0) = cube_mean.slice(0);
      c = p0ij_DR.slices(as<uvec>(s2));
      cube_mean = mean(c, 2);
      p0j_DR.row(1) = cube_mean.slice(0);
      c = p0ij_DR.slices(as<uvec>(s3));
      cube_mean = mean(c, 2);
      p0j_DR.row(2) = cube_mean.slice(0);
      c = p0ij_DR.slices(as<uvec>(s4));
      cube_mean = mean(c, 2);
      p0j_DR.row(3) = cube_mean.slice(0);
      c = p0ij_DR.slices(as<uvec>(s5));
      cube_mean = mean(c, 2);
      p0j_DR.row(4) = cube_mean.slice(0);
    }
    
    NumericVector E (n);
    E.fill(0);
    for (int i = 0; i<n; i++) {
      if(E1[i]>=E0[i]) {
        E[i] = 1;
      } else {
        E[i] = 0;
      }
    }
    NumericVector optimal (5);
    optimal.fill(0);
    NumericVector subE1 = E[s1];
    NumericVector subE2 = E[s2];
    NumericVector subE3 = E[s3];
    NumericVector subE4 = E[s4];
    NumericVector subE5 = E[s5];
    optimal[0] = mean(subE1);
    optimal[1] = mean(subE2);
    optimal[2] = mean(subE3);
    optimal[3] = mean(subE4);
    optimal[4] = mean(subE5);
    
    for (int g = 0; g<5; g++) {
      mat p1j = p1j_DR.row(g);
      mat p0j = p0j_DR.row(g);
      mat p (k+1,k+2);
      for (int m = 0; m<(k+1); m++) {
        p.row(m) = pj_R.row(g);
      }
      p1j_DR.row(g) = p1j/p;
      p0j_DR.row(g) = p0j/p;
    }
    
    arma::cube APCE = p1j_DR - p0j_DR;
    
    // save results
    out = List::create(Named("Pj_D1") = p1j_DR,
                            _["Pj_D0"] = p0j_DR,
                            _["APCE"] = APCE,
                            _["Pj_R"] = pj_R,
                            _["Optimal_Z"] = optimal,
                            _["Optimal_D"] = delta_star,
                            _["Optimal_D_ind"] = delta_star_i,
                            _["Utility_g_d_ind"] = g_d_i,
                            _["Utility_g_dmf_ind"] = g_dmf_i,
                            _["Pj_dmf"] = dmfj_DR);
  }
  
  
  return out;
}
