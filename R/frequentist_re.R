#' Test monotonicity with random effects
#' 
#' Test monotonicity using frequentist analysis with random effects for the hearing date of the case.
#' 
#' @param data A \code{data.frame} or \code{matrix} of which columns consists of pre-treatment covariates, a binary treatment (Z), an ordinal decision (D), and an outcome variable (Y). The column names of the latter three should be specified as "Z", "D", and "Y" respectively.
#' @param formula A formula of the model to fit.
#' 
#' @return Message indicating whether the monotonicity assumption holds.
#' 
#' @examples 
#' data(synth)
#' data(hearingdate_synth)
#' synth$CourtEvent_HearingDate = hearingdate_synth
#' TestMonotonicityRE(synth, formula = "Y ~ Sex + White + Age + 
#'                    CurrentViolentOffense + PendingChargeAtTimeOfOffense + 
#'                    PriorMisdemeanorConviction + PriorFelonyConviction + 
#'                    PriorViolentConviction + (1|CourtEvent_HearingDate) + D")
#' 
#' @importFrom lme4 glmer
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
TestMonotonicityRE = function(data,
                              formula){
  as.formula <- binomial <- NULL
  
  n = dim(data)[1]
  k=max(data$D)
  if (k<=1){ stop('Please use function for binary decision')}
  if(length(unique(data$Y))!=2){stop('Non-binary outcome')}
  
  data.factor=data
  data.factor$D=as.factor(data$D)
  formula = as.formula(formula)
  glm.e =glmer(formula = formula, family=binomial(link='probit'),data=data.factor)
  coefs = summary(glm.e)$coefficients
  delta = - coefs[grepl(paste0("(Intercept)|", paste0("D",1:k, collapse = "|")), rownames(coefs)),"Estimate"]
  Test.Pass = 1
  if (any(diff(delta)<0)){
    cat('Monotonicity Fails', delta)
    Test.Pass = 0}else{
      cat('Monotonicity Passes', delta)	
    }
  
}
#' Compute APCE using frequentist analysis with random effects
#' 
#' Estimate propensity score and use Hajek estimator to compute APCE. See S7 for more details.
#' 
#' @param data A \code{data.frame} or \code{matrix} of which columns consists of pre-treatment covariates, a binary treatment (Z), an ordinal decision (D), and an outcome variable (Y). The column names of the latter three should be specified as "Z", "D", and "Y" respectively.
#' @param formula A formula of the model to fit.
#' @param nAGQ Integer scalar - the number of points per axis for evaluating the adaptive Gauss-Hermite approximation to the log-likelihood. Defaults to 1, corresponding to the Laplace approximation. 
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
#' data(hearingdate_synth)
#' synth$CourtEvent_HearingDate = hearingdate_synth
#' freq_apce_re = CalAPCEipwRE(synth, formula = "Y ~ Sex + White + Age + 
#'                             CurrentViolentOffense + PendingChargeAtTimeOfOffense + 
#'                             PriorMisdemeanorConviction + PriorFelonyConviction + 
#'                             PriorViolentConviction + (1|CourtEvent_HearingDate) + D")
#' 
#' @importFrom lme4 glmer
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
CalAPCEipwRE <- function(data,
                         formula,
                         nAGQ = 1){
  as.formula <- binomial <- D <- predict <- NULL
  
  n = dim(data)[1]
  k=max(data$D)
  if (k<=1){ stop('Please use function for binary decision')}
  if(length(unique(data$Y))!=2){stop('Non-binary outcome')}
  ### (r+1)-th col of p.r: P(Y=1|D=r,X) for r = 0,...,k
  p.r=array(0,dim=c(n,k+1))
  
  data.factor=data
  data.factor$D=as.factor(data$D)
  formula = as.formula(formula)
  glm.e =glmer(formula = formula, family=binomial(link='probit'),data=data.factor, nAGQ = nAGQ)
  # if (any(diff(glm.e$coefficients[2:(k+1)])<0)){stop('The monotonicity is violated')}
  
  ### for comparison with MCMC algorithm
  coefs = summary(glm.e)$coefficients
  delta = - coefs[grepl(paste0("(Intercept)|", paste0("D",1:k, collapse = "|")), rownames(coefs)),"Estimate"]
  alpha = coefs[!grepl(paste0("(Intercept)|", paste0("D",1:k, collapse = "|")), rownames(coefs)),"Estimate"]
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
#' Bootstrap for estimating variance of APCE with random effects
#' 
#' Estimate variance of APCE for frequentist analysis with random effects using bootstrap. See S7 for more details.
#' 
#' @param data A \code{data.frame} or \code{matrix} of which columns consists of pre-treatment covariates, a binary treatment (Z), an ordinal decision (D), and an outcome variable (Y). The column names of the latter three should be specified as "Z", "D", and "Y" respectively.
#' @param rep Size of bootstrap
#' @param formula A formula of the model to fit.
#' @param CourtEvent_HearingDate The court event hearing date.
#' @param nAGQ Integer scalar - the number of points per axis for evaluating the adaptive Gauss-Hermite approximation to the log-likelihood. Defaults to 1, corresponding to the Laplace approximation. 
#' 
#' @return An object of class \code{list} with the following elements:
#' \item{P.D1.boot}{An array with dimension rep by (k+1) by (k+2) for quantity P(D(1)=d| R=r), dimension 1 is rep (size of bootstrap), dimension 2 is (k+1) values of D from 0 to k, dimension 3 is (k+2) values of R from 0 to k+1.}
#' \item{P.D0.boot}{An array with dimension rep by (k+1) by (k+2) for quantity P(D(0)=d| R=r).}
#' \item{APCE.boot}{An array with dimension rep by (k+1) by (k+2) for quantity P(D(1)=d| R=r)-P(D(0)=d| R=r).}
#' \item{P.R.boot}{An array with dimension rep by (k+2) for quantity P(R=r) for r from 0 to (k+1).}
#' 
#' @examples 
#' \dontrun{
#' data(synth)
#' data(hearingdate_synth)
#' synth$CourtEvent_HearingDate = hearingdate_synth
#' set.seed(123)
#' boot_apce_re = BootstrapAPCEipwRE(synth, rep = 10, formula = "Y ~ Sex + White + Age + 
#'                                   CurrentViolentOffense + PendingChargeAtTimeOfOffense + 
#'                                   PriorMisdemeanorConviction + PriorFelonyConviction + 
#'                                   PriorViolentConviction + (1|CourtEvent_HearingDate) + D", 
#'                                   CourtEvent_HearingDate = hearingdate_synth)
#' }
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr group_nest
#' @importFrom tidyr unnest
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
BootstrapAPCEipwRE <- function(data, rep=1000,
                               formula,
                               CourtEvent_HearingDate,
                               nAGQ = 1){
  n = dim(data)[1]
  p = dim(data)[2]-4
  k=max(data$D)
  if (k<=1){ stop('Please use function for binary decision')}
  
  #### bootstrap
  P.D1.boot=array(0,dim=c(rep,k+1,k+2))
  P.D0.boot=array(0,dim=c(rep,k+1,k+2))
  APCE.boot=array(0,dim=c(rep,k+1,k+2))
  P.R.boot = array(0,dim=c(rep,k+2))
  
  nest_data <- data %>%
    group_by(CourtEvent_HearingDate) %>%
    group_nest()
  m <- nrow(nest_data)
  
  for (i in 1:rep){
    if ( i %% 100 ==0) {print(paste('I am running the', i,'th repetion of boostrap'))}
    ind=sample(1:m,replace=TRUE)
    nest_data.boot=nest_data[ind,]
    data.boot=unnest(nest_data.boot, cols = "data")
    re.boot = CalAPCEipwRE(data.boot, formula, nAGQ)
    P.D1.boot[i,,]=re.boot$P.D1
    P.D0.boot[i,,]=re.boot$P.D0
    APCE.boot[i,,]=re.boot$APCE
    P.R.boot[i,]=re.boot$P.R
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
  
  return(list(P.D1.boot=P.D1.boot,P.D0.boot=P.D0.boot,APCE.boot=APCE.boot,P.R.boot=P.R.boot))
}
#' Bootstrap for estimating variance of APCE with random effects
#' 
#' Estimate variance of APCE for frequentist analysis with random effects using bootstrap. See S7 for more details.
#' 
#' @param data A \code{data.frame} or \code{matrix} of which columns consists of pre-treatment covariates, a binary treatment (Z), an ordinal decision (D), and an outcome variable (Y). The column names of the latter three should be specified as "Z", "D", and "Y" respectively.
#' @param rep Size of bootstrap
#' @param formula A formula of the model to fit.
#' @param nAGQ Integer scalar - the number of points per axis for evaluating the adaptive Gauss-Hermite approximation to the log-likelihood. Defaults to 1, corresponding to the Laplace approximation. 
#' @param size The number of parallel computing. The default is \code{5}.
#' 
#' @return An object of class \code{list} with the following elements:
#' \item{P.D1.boot}{An array with dimension rep by (k+1) by (k+2) for quantity P(D(1)=d| R=r), dimension 1 is rep (size of bootstrap), dimension 2 is (k+1) values of D from 0 to k, dimension 3 is (k+2) values of R from 0 to k+1.}
#' \item{P.D0.boot}{An array with dimension rep by (k+1) by (k+2) for quantity P(D(0)=d| R=r).}
#' \item{APCE.boot}{An array with dimension rep by (k+1) by (k+2) for quantity P(D(1)=d| R=r)-P(D(0)=d| R=r).}
#' \item{P.R.boot}{An array with dimension rep by (k+2) for quantity P(R=r) for r from 0 to (k+1).}
#' 
#' @examples 
#' \dontrun{
#' data(synth)
#' data(hearingdate_synth)
#' synth$CourtEvent_HearingDate = hearingdate_synth
#' set.seed(123)
#' boot_apce_re = BootstrapAPCEipwREparallel(synth, rep = 10, 
#'                                           formula = "Y ~ Sex + White + Age + 
#'                                           CurrentViolentOffense + PendingChargeAtTimeOfOffense + 
#'                                           PriorMisdemeanorConviction + PriorFelonyConviction + 
#'                                           PriorViolentConviction + (1|CourtEvent_HearingDate) + 
#'                                           D", size = 2)
#' }
#' 
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom purrr pluck
#' @importFrom abind abind
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
BootstrapAPCEipwREparallel <- function(data,
                                       rep=1000,
                                       formula,
                                       nAGQ = 1,
                                       size = 5) {
  j <- NULL
  numCores <- detectCores()
  registerDoParallel(numCores)
  lb = c(1, floor(rep/size)*1:(size-1) + 1)
  ub = c(floor(rep/size)*1:(size-1), rep)
  bootstrap_ls <- foreach (i=lb, j=ub) %dopar% {
    BootstrapAPCEipwRE(data,
                  rep=length(i:j),
                  formula,
                  nAGQ)
  }
  rs = list()
  l = names(bootstrap_ls[[1]])
  for (i in l) {
    rs[[i]]  = bootstrap_ls %>% 
      map(pluck(i)) %>%
      abind(along = 1)
  }
  return(rs)
}

