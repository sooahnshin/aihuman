#' Calculate the proportion of principal strata (R)
#' 
#' Calculate the proportion of each principal stratum (R).
#'
#' @param p.r.mcmc P.R.mcmc array generated from \code{CalAPCE} or \code{CalAPCEparallel}.
#' @param name.group A character vector including the labels of five subgroups.
#' 
#' @return A \code{data.frame} of the proportion of each principal stratum.
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#' 
#' @examples 
#' \dontrun{
#' data(synth)
#' sample_mcmc = AiEvalmcmc(data = synth, n.mcmc = 10)
#' subgroup_synth = list(1:nrow(synth),which(synth$Sex==0),which(synth$Sex==1),
#'                       which(synth$Sex==1&synth$White==0),which(synth$Sex==1&synth$White==1))
#' sample_apce = CalAPCEparallel(data = synth, mcmc.re = sample_mcmc, 
#'                               subgroup = subgroup_synth, size = 2)
#' CalPS(sample_apce[["P.R.mcmc"]])
#' }
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
CalPS <- function(p.r.mcmc, 
                  name.group = c("Overall", "Female", "Male", "Non-white\nMale", "White\nMale")) {
  subgroup <- value <- er <- NULL
  
  res <- rbind(apply(p.r.mcmc, c(2,3), mean),
               apply(p.r.mcmc, c(2,3), quantile, probs = 0.025),
               apply(p.r.mcmc, c(2,3), quantile, probs = 0.975))
  res <- as.data.frame(res)
  colnames(res) <- paste0("e",0:(ncol(res)-1))
  res <- res %>%
    mutate(subgroup = rep(name.group,3),
           stat = rep(c('mean','lb','ub'), each = 5)) %>%
    pivot_longer(!(subgroup|stat), names_to = "er", values_to = "value") %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    arrange(factor(subgroup, levels = name.group), er)
  
  return(res)
}

#' Plot the proportion of principal strata (R)
#' 
#' See Figure 3 for example.
#'
#' @param res A \code{data.frame} generated with \code{CalPS}.
#' @param y.min Minimum value of y-axis.
#' @param y.max Maximum value of y-axis.
#' @param col.values Color of point for each principal stratum.
#' @param label A logical argument whether to specify label of each principal stratum. The default is \code{TRUE}.
#' @param r.labels Label of each principal stratum.
#' @param label.position The position of labels.
#' @param top.margin Top margin of labels.
#' @param bottom.margin Bottom margin of labels.
#' @param label.size Size of label.
#' 
#' @return A ggplot.
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @import ggplot2
#' 
#' @examples 
#' \dontrun{
#' data(synth)
#' sample_mcmc = AiEvalmcmc(data = synth, n.mcmc = 10)
#' subgroup_synth = list(1:nrow(synth),which(synth$Sex==0),which(synth$Sex==1),
#'                       which(synth$Sex==1&synth$White==0),which(synth$Sex==1&synth$White==1))
#' sample_apce = CalAPCEparallel(data = synth, mcmc.re = sample_mcmc, 
#'                               subgroup = subgroup_synth, size = 2)
#' sample_ps = CalPS(sample_apce[["P.R.mcmc"]])
#' PlotPS(sample_ps, col.values = c("blue", "black", "red", "brown", "purple"), label = FALSE)
#' }
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
PlotPS <- function(res, 
                   y.min = 0, y.max = 0.75,
                   col.values = c("blue", "black", "red", "brown"),
                   label = TRUE,
                   r.labels = c("safe",
                               " easily             \n preventable    ",
                               "\n          preventable\n",
                               "  risky"),
                   label.position = c("top", "top", "top", "bottom"), 
                   top.margin = 0.02, bottom.margin = 0.02, 
                   label.size = 6.5) {
  idx <- er <- rgb <- lb <- ub <- NULL
  
  kp2 = length(col.values)
                   
  res <- res %>%
    mutate(idx = 1:nrow(res)+rep(0:4,each = kp2))
  
  p <- ggplot(res, aes(x=idx, y=mean, col = as.factor(er))) +
    geom_vline(xintercept = c(1:max(res$idx))[-res$idx], col = rgb(220,220,220, max = 255)) + 
    geom_point(shape = 18, size = 5) + theme_bw() + labs(x=NULL,y=NULL) +
    geom_errorbar(aes(ymin=lb, ymax=ub), width=.2) + 
    theme(legend.position="none", panel.grid = element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x = element_text(face="bold"),
          axis.text = element_text(size=20),
          legend.text = element_text(size=18),
          plot.title = element_text(face="bold", size = 25, hjust=0.5)) +
    ylim(y.min,y.max) +
    scale_color_manual(values=col.values) + 
    scale_x_continuous(breaks = c(0:4)*(kp2+1)+(kp2+1)/2, minor_breaks = NULL,
                       labels = unique(res$subgroup))
  
  if(label) {
    y.position = rep(0, kp2)
    for (i in 1:kp2) {
      if(label.position[i] == "top") {
        y.position[i] = res$ub[i]+top.margin
      } else {
        y.position[i] = res$lb[i]-bottom.margin
      }
      p <- p +
        annotate(geom="text", x=i, y=y.position[i], label=r.labels[i], size = label.size, col = col.values[i]) 
    }
  }
  
  return(p)
}
