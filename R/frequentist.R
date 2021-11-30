#' Test monotonicity 
#' 
#' Test monotonicity using frequentist analysis
#' 
#' @param data A \code{data.frame} or \code{matrix} of which columns consists of pre-treatment covariates, a binary treatment (Z), an ordinal decision (D), and an outcome variable (Y). The column names of the latter three should be specified as "Z", "D", and "Y" respectively.
#' 
#' @return Message indicating whether the monotonicity assumption holds.
#' 
#' @examples 
#' data(synth)
#' TestMonotonicity(synth)
#' 
#' @importFrom stats glm
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
TestMonotonicity= function(data){
  binomial <- NULL
  
  n = dim(data)[1]
  k=max(data$D)
  if (k<=1){ stop('Please use function for binary decision')}
  if(length(unique(data$Y))!=2){stop('Non-binary outcome')}
  
  data.factor=data
  data.factor$D=as.factor(data$D)
  glm.e =glm(Y~.-Z, family=binomial(link='probit'),data=data.factor)
  coefs = glm.e$coefficients
  delta = - coefs[grepl(paste0("(Intercept)|", paste0("D",1:k, collapse = "|")), names(coefs))]
  Test.Pass = 1
  if (any(diff(delta)<0)){
    cat('Monotonicity Fails', delta)
    Test.Pass = 0}else{
      cat('Monotonicity Passes', delta)	
    }
}
#' Compute APCE using frequentist analysis
#' 
#' Estimate propensity score and use Hajek estimator to compute APCE. See S7 for more details.
#' 
#' @param data A \code{data.frame} or \code{matrix} of which columns consists of pre-treatment covariates, a binary treatment (Z), an ordinal decision (D), and an outcome variable (Y). The column names of the latter three should be specified as "Z", "D", and "Y" respectively.
#' 
#' @return An object of class \code{list} with the following elements:
#' \item{P.D1}{An array with dimension (k+1) by (k+2) for quantity P(D(1)=d| R=r), dimension 1 is (k+1) values of D from 0 to k, dimension 2 is (k+2) values of R from 0 to k+1.}
#' \item{P.D0}{An array with dimension (k+1) by (k+2) for quantity P(D(0)=d| R=r).}
#' \item{APCE}{An array with dimension (k+1) by (k+2) for quantity P(D(1)=d| R=r)-P(D(0)=d| R=r).}
#' \item{P.R}{An array with dimension (k+2) for quantity P(R=r) for r from 0 to (k+1).}
#' \item{alpha}{An array with estimated alpha.}
#' \item{delta}{An array with estimated delta.}
#' 
#' @examples 
#' data(synth)
#' freq_apce = CalAPCEipw(synth)
#' 
#' @importFrom stats glm
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
CalAPCEipw= function(data){
  binomial <- D <- predict <- NULL
  
  n = dim(data)[1]
  k=max(data$D)
  if (k<=1){ stop('Please use function for binary decision')}
  if(length(unique(data$Y))!=2){stop('Non-binary outcome')}
  ### (r+1)-th col of p.r: P(Y=1|D=r,X) for r = 0,...,k
  p.r=array(0,dim=c(n,k+1))
  
  data.factor=data
  data.factor$D=as.factor(data$D)
  glm.e =glm(Y~.-Z, family=binomial(link='probit'),data=data.factor)
  # if (any(diff(glm.e$coefficients[2:(k+1)])<0)){stop('The monotonicity is violated')}
  
  ### for comparison with MCMC algorithm
  coefs = glm.e$coefficients
  delta = - coefs[grepl(paste0("(Intercept)|", paste0("D",1:k, collapse = "|")), names(coefs))]
  alpha = coefs[!grepl(paste0("(Intercept)|", paste0("D",1:k, collapse = "|")), names(coefs))]
  for (i in 0:k){
    data.new=data.frame(subset(data,select=-D),D=as.factor(rep(i,n)))
    p.r[,i+1]= predict(glm.e,newdata=data.new,type='response')
  }
  
  ### matrix for e_r(X)
  e.mat=array(0,dim=c(n,k+2))
  e.mat[,1]=1-p.r[,1]
  e.mat[,k+2]=p.r[,k+1]
  for (j in 1:k){
    e.mat[,j+1]=p.r[,j]-p.r[,j+1]
  }
  w.mat=array(0,dim=c(n,k+2))
  w.mat = t(t(e.mat)/ apply(e.mat,2,mean))
  P.D1 = array(0,dim=c(k+1,k+2))
  P.D0 = array(0,dim=c(k+1,k+2))
  APCE= array(0,dim=c(k+1,k+2))
  
  for ( j in 0:k){
    for (r in 0:(k+1)){
      P.D1[j+1,r+1] = mean(w.mat[data$Z==1,r+1]* (data$D[data$Z==1]==j))/mean(w.mat[data$Z==1,r+1])
      P.D0[j+1,r+1] = mean(w.mat[data$Z==0,r+1]* (data$D[data$Z==0]==j))/mean(w.mat[data$Z==0,r+1])
      APCE[j+1,r+1] = P.D1[j+1,r+1] -P.D0[j+1,r+1] 
    }
  }
  
  name.D=numeric(k+1)
  for (d in 0:k){
    name.D[d+1]=paste('D=',d,sep='')
  }
  name.R=numeric(k+2)
  for (r in 0:(k+1)){
    name.R[r+1]=paste('R=',r,sep='')
  }
  dimnames(P.D1)=list(name.D,name.R)  
  dimnames(P.D0)=list(name.D,name.R)    
  dimnames(APCE)=list(name.D,name.R)  
  
  e.mat = colMeans(e.mat)
  names(e.mat)=name.R
  
  return(list(P.D1=P.D1,P.D0=P.D0,APCE=APCE,P.R=e.mat,
              alpha=alpha,delta=delta))
  
}
#' Bootstrap for estimating variance of APCE
#' 
#' Estimate variance of APCE for frequentist analysis using bootstrap. See S7 for more details.
#' 
#' @param data A \code{data.frame} or \code{matrix} of which columns consists of pre-treatment covariates, a binary treatment (Z), an ordinal decision (D), and an outcome variable (Y). The column names of the latter three should be specified as "Z", "D", and "Y" respectively.
#' @param rep Size of bootstrap
#' 
#' @return An object of class \code{list} with the following elements:
#' \item{P.D1.boot}{An array with dimension rep by (k+1) by (k+2) for quantity P(D(1)=d| R=r), dimension 1 is rep (size of bootstrap), dimension 2 is (k+1) values of D from 0 to k, dimension 3 is (k+2) values of R from 0 to k+1.}
#' \item{P.D0.boot}{An array with dimension rep by (k+1) by (k+2) for quantity P(D(0)=d| R=r).}
#' \item{APCE.boot}{An array with dimension rep by (k+1) by (k+2) for quantity P(D(1)=d| R=r)-P(D(0)=d| R=r).}
#' \item{P.R.boot}{An array with dimension rep by (k+2) for quantity P(R=r) for r from 0 to (k+1).}
#' \item{alpha.boot}{An array with estimated alpha for each bootstrap.}
#' \item{delta.boot}{An array with estimated delta for each bootstrap.}
#' 
#' @examples 
#' data(synth)
#' set.seed(123)
#' boot_apce = BootstrapAPCEipw(synth, rep = 100)
#' 
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
BootstrapAPCEipw = function(data, rep=1000){
  n = dim(data)[1]
  p = dim(data)[2]-3
  k=max(data$D)
  if (k<=1){ stop('Please use function for binary decision')}
  
  #### bootstrap
  P.D1.boot=array(0,dim=c(rep,k+1,k+2))
  P.D0.boot=array(0,dim=c(rep,k+1,k+2))
  APCE.boot=array(0,dim=c(rep,k+1,k+2))
  P.R.boot = array(0,dim=c(rep,k+2))
  alpha.boot = array(0,dim=c(rep,p))
  delta.boot = array(0,dim=c(rep,k+1))
  for (i in 1:rep){
    if ( i %% 100 ==0) {print(paste('I am running the', i,'th repetion of boostrap'))}
    ind=sample(1:n,replace=TRUE)
    data.boot=data[ind,]
    re.boot = CalAPCEipw(data.boot)
    P.D1.boot[i,,]=re.boot$P.D1
    P.D0.boot[i,,]=re.boot$P.D0
    APCE.boot[i,,]=re.boot$APCE
    P.R.boot[i,]=re.boot$P.R
    alpha.boot[i,]=re.boot$alpha
    delta.boot[i,]=re.boot$delta
  }
  name.D=numeric(k+1)
  for (d in 0:k){
    name.D[d+1]=paste('D=',d,sep='')
  }
  name.R=numeric(k+2)
  for (r in 0:(k+1)){
    name.R[r+1]=paste('R=',r,sep='')
  }
  dimnames(P.D1.boot)=list(NULL,name.D,name.R)  
  dimnames(P.D0.boot)=list(NULL,name.D,name.R)    
  dimnames(APCE.boot)=list(NULL,name.D,name.R) 
  dimnames(P.R.boot)=list(NULL,name.R) 
  
  return(list(P.D1.boot=P.D1.boot,P.D0.boot=P.D0.boot,APCE.boot=APCE.boot,P.R.boot=P.R.boot,
              alpha.boot=alpha.boot,delta.boot=delta.boot))
}

#' Summary of APCE for frequentist analysis
#'
#' Summary of average principal causal effects (APCE) with ordinal decision with frequentist results.
#'
#' @param G1_est List generated from \code{CalAPCEipw} for the first subgroup.
#' @param G2_est List generated from \code{CalAPCEipw} for the second subgroup.
#' @param G3_est List generated from \code{CalAPCEipw} for the third subgroup.
#' @param G4_est List generated from \code{CalAPCEipw} for the fourth subgroup.
#' @param G5_est List generated from \code{CalAPCEipw} for the fifth subgroup.
#' @param G1_boot List generated from \code{BootstrapAPCEipw} for the first subgroup.
#' @param G2_boot List generated from \code{BootstrapAPCEipw} for the second subgroup.
#' @param G3_boot List generated from \code{BootstrapAPCEipw} for the third subgroup.
#' @param G4_boot List generated from \code{BootstrapAPCEipw} for the fourth subgroup.
#' @param G5_boot List generated from \code{BootstrapAPCEipw} for the fifth subgroup.
#' @param name.group A list of character vectors for the label of five subgroups.
#' 
#' @return A \code{data.frame} that consists of mean and quantiles (2.5%, 97.5%, 5%, 95%) of APCE (P(D(1)=d| R=r)-P(D(0)=d| R=r)) for each subgroup given specific value of D (decision) and R (principal strata).
#' 
#' @importFrom stats quantile
#' @importFrom dplyr bind_rows
#' 
#' @examples 
#' data(synth)
#' synth$SexWhite = synth$Sex * synth$White
#' freq_apce = CalAPCEipw(synth)
#' boot_apce = BootstrapAPCEipw(synth, rep = 10)
#' # subgroup analysis
#' data_s0 = subset(synth, synth$Sex==0,select=-c(Sex,SexWhite))
#' freq_s0 = CalAPCEipw(data_s0)
#' boot_s0 = BootstrapAPCEipw(data_s0, rep = 10)
#' data_s1 = subset(synth, synth$Sex==1,select=-c(Sex,SexWhite))
#' freq_s1 = CalAPCEipw(data_s1)
#' boot_s1 = BootstrapAPCEipw(data_s1, rep = 10)
#' data_s1w0 = subset(synth, synth$Sex==1&synth$White==0,select=-c(Sex,White,SexWhite))
#' freq_s1w0 = CalAPCEipw(data_s1w0)
#' boot_s1w0 = BootstrapAPCEipw(data_s1w0, rep = 10)
#' data_s1w1 = subset(synth, synth$Sex==1&synth$White==1,select=-c(Sex,White,SexWhite))
#' freq_s1w1 = CalAPCEipw(data_s1w1)
#' boot_s1w1 = BootstrapAPCEipw(data_s1w1, rep = 10)
#' 
#' freq_apce_summary <- APCEsummaryipw(freq_apce, freq_s0, freq_s1, freq_s1w0, freq_s1w1,
#'                                     boot_apce, boot_s0, boot_s1, boot_s1w0, boot_s1w0)
#' PlotAPCE(freq_apce_summary, y.max = 0.25, decision.labels = c("signature","small cash",
#'          "middle cash","large cash"), shape.values = c(16, 17, 15, 18), 
#'          col.values = c("blue", "black", "red", "brown", "purple"), label = FALSE)
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
APCEsummaryipw <- function(G1_est, G2_est, G3_est, G4_est, G5_est,
                    G1_boot, G2_boot, G3_boot, G4_boot, G5_boot,
                    name.group = c("Overall", "Female", "Male", "Non-white\nMale", "White\nMale")) {
  k=dim(G1_est$APCE)[1] - 1
  
  est_ls = list(G1_est, G2_est, G3_est, G4_est, G5_est)
  boot_ls = list(G1_boot, G2_boot, G3_boot, G4_boot, G5_boot)
  label = name.group
  
  res = list()
  for (i in 1:5) {
    est = est_ls[[i]]
    boot = boot_ls[[i]]
    ci = apply(boot$APCE.boot, c(2,3), quantile, probs = c(0.025, 0.975), na.rm = T)
    lb = c(ci[1,,])
    ub = c(ci[2,,])
    mean = c(est$APCE)
    res[[i]] = data.frame(subgroup = label[i],
                          D = rep(0:k, (k+2)),
                          R = rep(0:(k+1), each = (k+1)),
                          Mean = mean,
                          lb = lb,
                          ub = ub,
                          stringsAsFactors = F)
    colnames(res[[i]])[5:6] = c("2.5%", "97.5%")
  }
  res = bind_rows(res)
  
  return(res)
}
