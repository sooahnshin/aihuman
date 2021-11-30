#' Calculate the delta given the principal stratum
#' 
#' Calculate the maximal deviation of decisions probability among the distributions for different groups (delta) given the principal stratum (R).
#'
#' @param r The given principal stratum.
#' @param k The maximum decision (e.g., largest bail amount).
#' @param pd0 P.D0.mcmc array generated from \code{CalAPCE} or \code{CalAPCEparallel}.
#' @param pd1 P.D1.mcmc array generated from \code{CalAPCE} or \code{CalAPCEparallel}.
#' @param attr The index of subgroups (within the output of CalAPCE/CalAPCEparallel) that corresponds to the protected attributes.
#' 
#' @return A \code{data.frame} of the delta.
#' 
#' @examples 
#' \dontrun{
#' data(synth)
#' subgroup_synth = list(1:nrow(synth), which(synth$Sex==0), which(synth$Sex==1), 
#'                       which(synth$Sex==1&synth$White==0), which(synth$Sex==1&synth$White==1))
#' sample_mcmc = AiEvalmcmc(data = synth, n.mcmc = 10)
#' sample_apce = CalAPCEparallel(data = synth, mcmc.re = sample_mcmc, subgroup = subgroup_synth, 
#'                               burnin = 0, size = 2)
#' CalDelta(0, 3, sample_apce[["P.D0.mcmc"]], sample_apce[["P.D1.mcmc"]], c(2,3))
#' }
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
CalDelta <- function(r, k, pd0, pd1, attr) {
  p0 = pd0[,attr,,r+1]
  p1 = pd1[,attr,,r+1]
  for (j in 1:(k+1)){
    p0[,,j] = apply(p0[,,j:(k+1)],c(1,2),sum)
    p1[,,j] = apply(p1[,,j:(k+1)],c(1,2),sum)
  }
  cand0 = apply(p0,c(2,3),mean)
  cand1 = apply(p1,c(2,3),mean)
  
  a0 = a1 = rep(NA, k+1)
  for(i in 1:(k+1)) {
    s1 = which.max(cand0[,i])
    s2 = which.min(cand0[,i])
    a0[i] = cand0[s1,i]-cand0[s2,i]
    
    s1 = which.max(cand1[,i])
    s2 = which.min(cand1[,i])
    a1[i] = cand1[s1,i]-cand1[s2,i]
  }
  d0 = which.max(a0)
  d1 = which.max(a1)
  
  a01 = which.max(cand0[,d0])
  a02 = which.min(cand0[,d0])
  a11 = which.max(cand1[,d1])
  a12 = which.min(cand1[,d1])
  
  delta0 = abs(p0[,a01,d0]-p0[,a02,d0]) 
  delta1 = abs(p1[,a11,d1]-p1[,a12,d1]) 
  
  out = data.frame(mean = apply(cbind(delta0, delta1, delta1-delta0), 2, mean), 
                   lb = apply(cbind(delta0, delta1, delta1-delta0), 2, quantile, probs = 0.025),
                   ub = apply(cbind(delta0, delta1, delta1-delta0), 2, quantile, probs = 0.975), 
                   s1 = c(names(which.max(cand0[,d0])), names(which.max(cand1[,d1])), NA), 
                   s2 = c(names(which.min(cand0[,d0])), names(which.min(cand1[,d1])), NA),
                   d = c(d0-1, d1-1, NA),
                   param = c("Delta0", "Delta1", "Delta"),
                   R = r)
  
  return(out)
}

#' Calculate the principal fairness
#' 
#' See Section 3.6 for more details.
#'
#' @param apce The list generated from \code{CalAPCE} or \code{CalAPCEparallel}.
#' @param attr The index of subgroups (within the output of CalAPCE/CalAPCEparallel) that corresponds to the protected attributes.
#' 
#' @return A \code{data.frame} of the delta.
#' 
#' @examples 
#' \dontrun{
#' data(synth)
#' subgroup_synth = list(1:nrow(synth), which(synth$Sex==0), which(synth$Sex==1), 
#'                       which(synth$Sex==1&synth$White==0), which(synth$Sex==1&synth$White==1))
#' sample_mcmc = AiEvalmcmc(data = synth, n.mcmc = 10)
#' sample_apce = CalAPCEparallel(data = synth, mcmc.re = sample_mcmc, subgroup = subgroup_synth, 
#'                               burnin = 0, size = 2)
#' CalFairness(sample_apce)
#' }
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
CalFairness <- function(apce, attr = c(2,3)) {
  pd1 = apce$P.D1.mcmc # P(D(1)=d| R=r)
  pd0 = apce$P.D0.mcmc # P(D(0)=d| R=r)
  
  k = dim(apce$APCE.mcmc)[3]-1
  
  res = list()
  for (r in 0:(k+1)) {
    res[[r+1]] = CalDelta(r, k, pd0, pd1, attr)
  }
  
  res = do.call(rbind, res)
  
  return(res)
}

#' Plot the principal fairness
#' 
#' See Figure 5 for example.
#'
#' @param res The data frame generated from \code{CalFairness}.
#' @param top.margin The index of subgroups (within the output of CalAPCE/CalAPCEparallel) that corresponds to the protected attributes.
#' @param y.max Maximum value of y-axis.
#' @param y.min Minimum value of y-axis.
#' @param r.labels Label of each principal stratum.
#' @param label A logical argument whether to specify label.
#' 
#' @return A \code{data.frame} of the delta.
#' 
#' @examples 
#' \dontrun{
#' data(synth)
#' subgroup_synth = list(1:nrow(synth), which(synth$Sex==0), which(synth$Sex==1), 
#'                       which(synth$Sex==1&synth$White==0), which(synth$Sex==1&synth$White==1))
#' sample_mcmc = AiEvalmcmc(data = synth, n.mcmc = 10)
#' sample_apce = CalAPCEparallel(data = synth, mcmc.re = sample_mcmc, subgroup = subgroup_synth, 
#'                               burnin = 0, size = 2)
#' sample_fair = CalFairness(sample_apce)
#' PlotFairness(sample_fair, y.max = 0.4, y.min = -0.4, r.labels = c("Safe", "Preventable 1", 
#'              "Preventable 2", "Preventable 3", "Risky"))
#' }
#' 
#' @import ggplot2
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
PlotFairness <- function(res, 
                         top.margin = 0.01, 
                         y.max = 0.2, y.min = -0.1,
                         r.labels = c("Safe",
                                      "Easily\nPreventable",
                                      "Preventable",
                                      "Risky"),
                         label = TRUE) {
  idx <- param <- rgb <- lb <- ub <- NULL
  
  kp1 = max(res$R)
  res$idx = 1:nrow(res)+rep(0:kp1,each = 3)
  
  p <- ggplot(res, aes(x=idx, y=mean, shape = as.factor(param))) +
    geom_vline(xintercept = c(1:max(res$idx))[-res$idx], col = rgb(220,220,220, max = 255)) + 
    geom_hline(yintercept = 0,  col = "darkgrey", lty = "dashed") +
    geom_point(size = 5) + theme_bw() + labs(x=NULL,y=NULL) +
    geom_errorbar(aes(ymin=lb, ymax=ub), width=.2) + 
    theme(legend.position="none", panel.grid = element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x = element_text(face="bold"),
          axis.text = element_text(size=20),
          legend.text = element_text(size=18),
          plot.title = element_text(face="bold", size = 25, hjust=0.5),
          axis.title=element_text(size=20)) +
    ylim(y.min,y.max) +
    scale_x_continuous(breaks = (0:kp1)*4+2, minor_breaks = NULL,
                       labels = r.labels) 
  
  if (label) {
    y.position = rep(0, 3)
    for (i in 1:3) {
      y.position[i] = res$ub[i]+top.margin
    }
    p <- p + 
      annotate(geom="text", x=1, y=y.position[1], label="Delta(0)", size = 5, parse = T) + 
      annotate(geom="text", x=2, y=y.position[2], label="Delta(1)", size = 5, parse = T) + 
      annotate(geom="text", x=3, y=y.position[3], label="Delta(1)-Delta(0)", size = 5, parse = T)
  }
  
  return(p)
}

