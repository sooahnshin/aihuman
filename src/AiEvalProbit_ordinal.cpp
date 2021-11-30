#include <RcppArmadillo.h>
#define NDEBUG
#include <RcppEigen.h>

using namespace Rcpp;
using namespace arma;
using namespace Eigen;

mat rmvnorm(int n, const vec& mu, const mat& sigma) {
  int ncols = sigma.n_cols;
  mat Y = randn(n, ncols);
  return repmat(mu, 1, n).t() + Y * chol(sigma);
}

MatrixXd matmult(NumericVector a){
  Map<VectorXd> A = as<Map<VectorXd> >(a);
  MatrixXd B = A * A.transpose();
  return B;
}

double rtruncnorm(double lb, double ub, double mu, double sigma) {
  double z;
  double stlb = (lb-mu)/sigma;  
  double stub = (ub-mu)/sigma; 
  if(stlb >= stub) {
    stop("TruncNorm: lower bound is greater than upper bound\n");
  }
  double tol=2.0;
  double temp;
  double M;
  double u;
  double exp_par;
  int flag=0;  
  if(stub<=-tol){
    flag=1;
    temp=stub;
    stub=-stlb;
    stlb=-temp;
  }
  if(stlb>=tol){
    exp_par=stlb;
    while(R::pexp(stub,1/exp_par,1,0) - R::pexp(stlb,1/exp_par,1,0) < 0.000001) {
      exp_par/=2.0;
    }
    if(R::dnorm(stlb,0,1,1) - R::dexp(stlb,1/exp_par,1) >= R::dnorm(stub,0,1,1) - R::dexp(stub,1/exp_par,1)) {
      M=exp(R::dnorm(stlb,0,1,1) - R::dexp(stlb,1/exp_par,1));
    } else {
      M=exp(R::dnorm(stub,0,1,1) - R::dexp(stub,1/exp_par,1));
    }
    u=R::unif_rand();
    z=-log(1-u*(R::pexp(stub,1/exp_par,1,0)-R::pexp(stlb,1/exp_par,1,0))-R::pexp(stlb,1/exp_par,1,0))/exp_par;
    while(unif_rand() > exp(R::dnorm(z,0,1,1)-R::dexp(z,1/exp_par,1))/M) {
      u=R::unif_rand();
      z=-log(1-u*(R::pexp(stub,1/exp_par,1,0)-R::pexp(stlb,1/exp_par,1,0))-R::pexp(stlb,1/exp_par,1,0))/exp_par;
    }
    if(flag==1) z=-z;
  } else {
    z=R::norm_rand();
    while((z<stlb)|(z>stub)) {
      z=R::norm_rand();
    }
  }
  return z*sigma + mu;
}

// [[Rcpp::export(.AiEvalProbitOrdinal)]]
NumericMatrix AiEvalProbitOrdinal(NumericVector Z,
                                  NumericVector D,
                                  NumericVector Y,
                                  NumericMatrix X,
                                  NumericMatrix D_factor,
                                  NumericMatrix ZX,
                                  double rho,
                                  Nullable<NumericMatrix> RSigma0_beta_inv,
                                  Nullable<NumericMatrix> RSigma0_alpha_inv,
                                  Nullable<double> Rsigma0,
                                  Nullable<NumericVector> Rbeta,
                                  Nullable<NumericVector> Ralpha,
                                  Nullable<NumericVector> Rtheta,
                                  Nullable<NumericVector> Rdelta,
                                  int n_mcmc, 
                                  bool verbose, 
                                  int out_length,
                                  bool Zbeta_off,
                                  bool theta01_off,
                                  int lZX,
                                  arma::mat C){
  // default parameters
  NumericVector dim = X.attr("dim");
  int p = dim[1];
  int n = dim[0];
  int k = unique(as<NumericVector>(D)).length()-1;
  int tk = unique(as<NumericVector>(D)).length()-1;
  if (!theta01_off) {
    tk = 2*tk;
  }
  if(k<2){
    stop("Please use functions for binary decision\n");
  }
  mat Sigma0_beta_inv = diagmat(ones(p+1+lZX)*0.01);
  mat Sigma0_alpha_inv = diagmat(ones(p)*0.01);
  double sigma0 = 10;
  NumericVector alpha (p);
  NumericVector beta (p+1+lZX);
  NumericVector theta (2*k);
  theta[0] = -1;
  if (!theta01_off) {
    theta[k] = -1;
    for (int i = 1; i < k; i++) {
      theta[i] = rtruncnorm(theta[i-1], R_PosInf, 0, 1);
      theta[k+i] = theta[i];
    }
  } else {
    for (int i = 1; i < k; i++) {
      theta[i] = rtruncnorm(theta[i-1], R_PosInf, 0, 1);
    }
  }
  
  NumericVector delta (k+1);
  delta[0] = -1;
  for (int i = 1; i < (k+1); i++) {
    delta[i] = rtruncnorm(delta[i-1], R_PosInf, 0, 1);
  }
  // generating the prior parameters
  if(!RSigma0_beta_inv.isNull()) {
    Sigma0_beta_inv = as<mat>(RSigma0_beta_inv);
  }
  if(!RSigma0_alpha_inv.isNull()) {
    Sigma0_alpha_inv = as<mat>(RSigma0_alpha_inv);
  }
  if(!Rsigma0.isNull()) {
    sigma0 = as<double>(Rsigma0);
  }
  // generating initial values
  if(!Rbeta.isNull()) {
    beta = as<NumericVector>(Rbeta);
  }
  if(!Ralpha.isNull()) {
    alpha = as<NumericVector>(Ralpha);
  }
  if(!Rtheta.isNull()) {
    theta = as<NumericVector>(Rtheta);
  }
  if(!Rdelta.isNull()) {
    delta = as<NumericVector>(Rdelta);
  }
  // transform data
  // construct W
  NumericMatrix W (n,k+1);
  int Wsize = W.nrow() * W.ncol();
  for (int i = 0; i < Wsize; i++) {
    W[i] = 1;
  }
  for (int i = 0; i < n; i++) {
    NumericVector D_cumsum = cumsum(D_factor(i,_));
    W(i,_) = W(i,_)-D_cumsum+D_factor(i,_);
  }
  NumericVector D_star (n);
  NumericVector Y_star (n);
  for (int i = 0; i < n; i++) {
    if ((D[i]>0) & (D[i]<k)) {
      D_star[i] = rtruncnorm(theta[D[i]-1], theta[D[i]], 0, 1);
    } else {
      if (D[i]==0) {
        D_star[i] = rtruncnorm(R_NegInf, theta[0], 0, 1);
      }
      if (D[i]==k) {
        D_star[i] = rtruncnorm(theta[k-1], R_PosInf, 0, 1);
      }
    }
  }
  // the mcmc results
  NumericMatrix BETA_mcmc(n_mcmc+1, p+1+lZX);
  NumericMatrix ALPHA_mcmc(n_mcmc+1, p);
  NumericMatrix THETA_mcmc(n_mcmc+1, tk);
  NumericMatrix DELTA_mcmc(n_mcmc+1, k+1);
  // save initial values
  BETA_mcmc(0,_) = beta;
  ALPHA_mcmc(0,_) = alpha;
  THETA_mcmc(0,_) = theta;
  DELTA_mcmc(0,_) = delta;
  NumericVector ninf = NumericVector::create(R_NegInf, 0);
  NumericVector pinf = NumericVector::create(0, R_PosInf);
  List ZXXZ(n);
  List XX(n);
  List WW(n);
  for (int i = 0; i<n; i++) {
    ZXXZ[i] = matmult(ZX(i,_));
    XX[i] = matmult(X(i,_));
    WW[i] = matmult(W(i,_));
  }
  // Gibbs Sampler
  if (!theta01_off) {
    for (int mcmc = 1; mcmc <= n_mcmc; mcmc++) {
      // imputation step
      NumericVector theta1_extend = theta[Range(0,k-1)];
      NumericVector theta0_extend = theta[Range(k,2*k-1)];
      theta1_extend.push_front(R_NegInf);
      theta1_extend.push_back(R_PosInf);
      theta0_extend.push_front(R_NegInf);
      theta0_extend.push_back(R_PosInf);
      for (int i=0; i<n; i++) {
        // sample D_star
        double mu_Ystari = - sum(D_factor(i,_)*delta) + sum(X(i,_)*alpha);
        double mu_Dstari = sum(ZX(i,_)*beta);
        double y_lower = ninf[Y[i]];
        double y_upper = pinf[Y[i]];
        Y_star[i] = rtruncnorm(y_lower, y_upper, mu_Ystari-rho*(D_star[i]-mu_Dstari), sqrt(1-pow(rho, 2.0)));
        if (Z[i]==1) {
          D_star[i] = rtruncnorm(theta1_extend[D[i]], theta1_extend[D[i]+1], mu_Dstari-rho*(Y_star[i]-mu_Ystari), sqrt(1-pow(rho, 2.0)));
        } else {
          D_star[i] = rtruncnorm(theta0_extend[D[i]], theta0_extend[D[i]+1], mu_Dstari-rho*(Y_star[i]-mu_Ystari), sqrt(1-pow(rho, 2.0)));
        }
      }
      // Posterior step
      // sample beta
      mat Sigma_beta = zeros(p+1+lZX,p+1+lZX);
      vec mu_beta(p+1+lZX);
      mu_beta.zeros();
      for (int i = 0; i<n; i++) {
        mat M = ZXXZ[i];
        Sigma_beta += M;
        mu_beta += ZX(i,_)*(D_star[i]-rho*(Y_star[i]+sum(D_factor(i,_)*delta)-sum(X(i,_)*alpha)));
      }
      mu_beta /= (1-pow(rho, 2.0));
      Sigma_beta = Sigma_beta/(1-pow(rho, 2.0)) + Sigma0_beta_inv;
      Sigma_beta = inv(Sigma_beta);
      mu_beta = Sigma_beta*mu_beta;
      beta = rmvnorm(1, mu_beta, Sigma_beta);
      BETA_mcmc(mcmc,_) = beta;
      // sample alpha
      mat Sigma_alpha = zeros(p,p);
      vec mu_alpha(p);
      mu_alpha.zeros();
      for (int i = 0; i<n; i++) {
        mat M = XX[i];
        Sigma_alpha += M;
        mu_alpha += X(i,_)*(-rho*(D_star[i]-sum(ZX(i,_)*beta))+Y_star[i]+sum(D_factor(i,_)*delta));
      }
      mu_alpha /= (1-pow(rho, 2.0));
      Sigma_alpha = Sigma_alpha/(1-pow(rho, 2.0)) + Sigma0_alpha_inv;
      Sigma_alpha = inv(Sigma_alpha);
      mu_alpha = Sigma_alpha*mu_alpha;
      alpha = rmvnorm(1, mu_alpha, Sigma_alpha);
      ALPHA_mcmc(mcmc,_) = alpha;
      // sample delta
      mat Sigma_delta = zeros(k+1,k+1);
      vec mu_delta(k+1);
      mu_delta.zeros();
      for (int i = 0; i<n; i++) {
        mat M = WW[i];
        Sigma_delta += M;
        mu_delta += W(i,_)*(-Y_star[i]+sum(X(i,_)*alpha)+rho*(D_star[i]- sum(ZX(i,_)*beta)));
      }
      mu_delta /= (1-pow(rho, 2.0));
      Sigma_delta = Sigma_delta/(1-pow(rho, 2.0)) + C/pow(sigma0, 2.0); // C^\transpose * C/sigma_0^2 
      Sigma_delta = inv(Sigma_delta);
      mu_delta = Sigma_delta*mu_delta;
      for (int j = 0; j<(k+1); j++) {
        colvec delta_rj(delta.begin(),delta.size(),false);
        delta_rj.shed_row(j);
        colvec mu_delta_rj = mu_delta;
        mu_delta_rj.shed_row(j);
        rowvec Sig_delta_cj = Sigma_delta.row(j);
        Sig_delta_cj.shed_col(j);
        colvec Sig_delta_rj = Sigma_delta.col(j);
        Sig_delta_rj.shed_row(j);
        mat S = Sigma_delta;
        S.shed_col(j);
        S.shed_row(j);
        S = inv(S);
        double mu_deltaj = mu_delta[j] + sum(Sig_delta_cj*S*(delta_rj-mu_delta_rj));
        double sigma2_deltaj = Sigma_delta(j,j) - sum(Sig_delta_cj*S*Sig_delta_rj);
        if (j==0) {
          delta[j] = rnorm(1, mu_deltaj, sqrt(sigma2_deltaj))[0];
        } else {
          delta[j] = rtruncnorm(0,R_PosInf,mu_deltaj,sqrt(sigma2_deltaj));
        }
      }
      NumericVector delta_cumsum = cumsum(delta);
      DELTA_mcmc(mcmc,_) = delta_cumsum;
      // sample theta
      for (int j=0; j<k; j++) {
        NumericVector theta1_lb = NumericVector::create(theta1_extend[j]);
        NumericVector theta1_ub = NumericVector::create(theta1_extend[j+2]);
        NumericVector theta0_lb = NumericVector::create(theta0_extend[j]);
        NumericVector theta0_ub = NumericVector::create(theta0_extend[j+2]);
        for (int i=0; i<n; i++) {
          if (Z[i]==1) {
            if (D[i]==j){
              theta1_lb.push_back(D_star[i]);
            } else if (D[i]==(j+1)){
              theta1_ub.push_back(D_star[i]);
            }
          } else {
            if (D[i]==j){
              theta0_lb.push_back(D_star[i]);
            } else if (D[i]==(j+1)){
              theta0_ub.push_back(D_star[i]);
            }
          }
        }
        theta[j] = rtruncnorm(max(theta1_lb),min(theta1_ub),0,sigma0);
        theta[k+j] = rtruncnorm(max(theta0_lb),min(theta0_ub),0,sigma0);
      }
      THETA_mcmc(mcmc,_) = theta;
      // print the mcmc process
      if (mcmc % out_length == 0) {
        Rprintf("%i/%i done.\n", mcmc, n_mcmc);
      }
    }
  } else {
    for (int mcmc = 1; mcmc <= n_mcmc; mcmc++) {
      // imputation step
      NumericVector theta_extend = theta[Range(0,k-1)];
      theta_extend.push_front(R_NegInf);
      theta_extend.push_back(R_PosInf);
      for (int i=0; i<n; i++) {
        // sample D_star
        double mu_Ystari = - sum(D_factor(i,_)*delta) + sum(X(i,_)*alpha);
        double mu_Dstari = sum(ZX(i,_)*beta);
        double y_lower = ninf[Y[i]];
        double y_upper = pinf[Y[i]];
        Y_star[i] = rtruncnorm(y_lower, y_upper, mu_Ystari-rho*(D_star[i]-mu_Dstari), sqrt(1-pow(rho, 2.0)));
        D_star[i] = rtruncnorm(theta_extend[D[i]], theta_extend[D[i]+1], mu_Dstari-rho*(Y_star[i]-mu_Ystari), sqrt(1-pow(rho, 2.0)));
      }
      // Posterior step
      // sample beta
      mat Sigma_beta = zeros(p+1+lZX,p+1+lZX);
      vec mu_beta(p+1+lZX);
      mu_beta.zeros();
      for (int i = 0; i<n; i++) {
        mat M = ZXXZ[i];
        Sigma_beta += M;
        mu_beta += ZX(i,_)*(D_star[i]-rho*(Y_star[i]+sum(D_factor(i,_)*delta)-sum(X(i,_)*alpha)));
      }
      mu_beta /= (1-pow(rho, 2.0));
      Sigma_beta = Sigma_beta/(1-pow(rho, 2.0)) + Sigma0_beta_inv;
      Sigma_beta = inv(Sigma_beta);
      mu_beta = Sigma_beta*mu_beta;
      beta = rmvnorm(1, mu_beta, Sigma_beta);
      BETA_mcmc(mcmc,_) = beta;
      // sample alpha
      mat Sigma_alpha = zeros(p,p);
      vec mu_alpha(p);
      mu_alpha.zeros();
      for (int i = 0; i<n; i++) {
        mat M = XX[i];
        Sigma_alpha += M;
        mu_alpha += X(i,_)*(-rho*(D_star[i]- sum(ZX(i,_)*beta))+Y_star[i]+sum(D_factor(i,_)*delta));
      }
      mu_alpha /= (1-pow(rho, 2.0));
      Sigma_alpha = Sigma_alpha/(1-pow(rho, 2.0)) + Sigma0_alpha_inv;
      Sigma_alpha = inv(Sigma_alpha);
      mu_alpha = Sigma_alpha*mu_alpha;
      alpha = rmvnorm(1, mu_alpha, Sigma_alpha);
      ALPHA_mcmc(mcmc,_) = alpha;
      // sample delta
      mat Sigma_delta = zeros(k+1,k+1);
      vec mu_delta(k+1);
      mu_delta.zeros();
      for (int i = 0; i<n; i++) {
        mat M = WW[i];
        Sigma_delta += M;
        mu_delta += W(i,_)*(-Y_star[i]+sum(X(i,_)*alpha)+rho*(D_star[i]- sum(ZX(i,_)*beta)));
      }
      mu_delta /= (1-pow(rho, 2.0));
      Sigma_delta = Sigma_delta/(1-pow(rho, 2.0))  + C/pow(sigma0, 2.0); // C^\transpose * C/sigma_0^2 
      Sigma_delta = inv(Sigma_delta);
      mu_delta = Sigma_delta*mu_delta;
      for (int j = 0; j<(k+1); j++) {
        colvec delta_rj(delta.begin(),delta.size(),false);
        delta_rj.shed_row(j);
        colvec mu_delta_rj = mu_delta;
        mu_delta_rj.shed_row(j);
        rowvec Sig_delta_cj = Sigma_delta.row(j);
        Sig_delta_cj.shed_col(j);
        colvec Sig_delta_rj = Sigma_delta.col(j);
        Sig_delta_rj.shed_row(j);
        mat S = Sigma_delta;
        S.shed_col(j);
        S.shed_row(j);
        S = inv(S);
        double mu_deltaj = mu_delta[j] + sum(Sig_delta_cj*S*(delta_rj-mu_delta_rj));
        double sigma2_deltaj = Sigma_delta(j,j) - sum(Sig_delta_cj*S*Sig_delta_rj);
        if (j==0) {
          delta[j] = rnorm(1, mu_deltaj, sqrt(sigma2_deltaj))[0];
        } else {
          delta[j] = rtruncnorm(0,R_PosInf,mu_deltaj,sqrt(sigma2_deltaj));
        }
      }
      NumericVector delta_cumsum = cumsum(delta);
      DELTA_mcmc(mcmc,_) = delta_cumsum;
      // sample theta
      for (int j=0; j<k; j++) {
        NumericVector theta_lb = NumericVector::create(theta_extend[j]);
        NumericVector theta_ub = NumericVector::create(theta_extend[j+2]);
        for (int i=0; i<n; i++) {
          if (D[i]==j){
            theta_lb.push_back(D_star[i]);
          } else if (D[i]==(j+1)){
            theta_ub.push_back(D_star[i]);
          }
        }
        theta[j] = rtruncnorm(max(theta_lb),min(theta_ub),0,sigma0);
      }
      THETA_mcmc(mcmc,_) = theta;
      // print the mcmc process
      if (mcmc % out_length == 0) {
        Rprintf("%i/%i done.\n", mcmc, n_mcmc);
      }
    }
  }
  // save results
  NumericMatrix out(n_mcmc+1, (p+p+k+tk+2+lZX));
  for (int j = 0; j < (p+p+k+tk+2+lZX); j++) {
    if (j < p+1+lZX) {
      out(_, j) = BETA_mcmc(_, j);
    } else if (j < p+1+lZX+tk){
      out(_, j) = THETA_mcmc(_, j - (p+1+lZX));
    } else if (j < p+1+lZX+tk+p){
      out(_, j) = ALPHA_mcmc(_, j - (p+1+lZX+tk));
    } else {
      out(_, j) = DELTA_mcmc(_, j - (p+1+lZX+tk+p));
    }
  }
  return out;
}
