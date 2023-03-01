#' Plot utility difference
#'
#' See Figure 7 for example.
#' 
#' @param res The data frame generated from \code{CalUtilityDiff}.
#' @param idx The row index of observations to be included. The default is all the observations from the data.
#' 
#' @return A ggplot.
#' 
#' @examples 
#' \donttest{
#' data(synth)
#' sample_mcmc = AiEvalmcmc(data = synth, n.mcmc = 10)
#' synth_dmf = sample(0:1, nrow(synth), replace = TRUE) # random dmf recommendation
#' sample_utility = CalOptimalDecision(data = synth, mcmc.re = sample_mcmc, 
#'                                     c0.ls = seq(0,5,1), c1.ls = seq(0,5,1), 
#'                                     dmf = synth_dmf, size = 1) # adjust the size
#' PlotUtilityDiff(sample_utility)
#' }
#' 
#' @import ggplot2
#' @importFrom metR geom_text_contour
#' @importFrom metR label_placer_flattest
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
PlotUtilityDiff <- function(res,
                            idx = NULL) {
  c0 <- c1 <- diff_utility <- NULL
  
  if(!is.null(idx)) {
    res = res[res$idx %in% idx,]
  }
  
  dat <- res %>%
    group_by(c0, c1) %>%
    summarise(diff_utility = mean(diff_utility)) 
  
  p <- ggplot(res, aes(c0, c1, z=diff_utility, fill= diff_utility)) + 
    geom_tile(alpha = 0.8)+
    scale_fill_distiller(limits = c(0,1),
                         palette = "Greys") +
    stat_contour(color = "black", alpha = 0.8)+
    geom_text_contour(color = "black", label.placer = label_placer_flattest()) +
    theme_bw() +
    theme(axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), panel.border = element_blank(),
          legend.position = "none",
          panel.background = element_blank(), 
          axis.text = element_text(size=18),
          legend.title = element_text(size=12),
          legend.text = element_text(size=12),
          plot.title = element_text(face="bold", size = 20, hjust=0.5),
          plot.subtitle = element_text(face="bold", size = 15, hjust=0.5),
          axis.title=element_text(size=18)) + 
    labs(x = expression("Cost of outcome ("*c[0]*")"),
         y = expression("Cost of unnecessarily harsh decision ("*c[1]*")"))
  
  return(p)
}
#' Plot utility difference with 95\% confidence interval
#'
#' See Figure S17 for example.
#' 
#' @param res The second data frame (res.mcmc) generated from \code{CalUtilityDiff(include.utility.diff.mcmc = TRUE)}.
#' 
#' @return A ggplot.
#' 
#' @examples 
#' \donttest{
#' data(synth)
#' sample_mcmc = AiEvalmcmc(data = synth, n.mcmc = 10)
#' synth_dmf = sample(0:1, nrow(synth), replace = TRUE) # random dmf recommendation
#' sample_utility = CalOptimalDecision(data = synth, mcmc.re = sample_mcmc, 
#'                                     c0.ls = seq(0,5,1), c1.ls = seq(0,5,1), 
#'                                     dmf = synth_dmf, size = 1, # adjust the size
#'                                     include.utility.diff.mcmc = TRUE)
#' PlotUtilityDiffCI(sample_utility$res.mcmc)
#' }
#' 
#' @import ggplot2
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
PlotUtilityDiffCI <- function(res) {
  c0 <- mean_treated_utility_diff <- atop <- lb_treated_utility_diff <- ub_treated_utility_diff <- mean_control_utility_diff <- lb_control_utility_diff <- ub_control_utility_diff <- NULL
  
  p1 <- ggplot(res, aes(x=c0, y=mean_treated_utility_diff)) +
    geom_hline(yintercept = 0,  col = "darkgrey", lty = "dashed") + 
    geom_point()+facet_wrap(~ c1, labeller = label_bquote(atop("Cost of unnecessarily","harsh decision ("*c[1]*")"~"="~.(c1)))) +
    geom_errorbar(aes(ymin=lb_treated_utility_diff, ymax=ub_treated_utility_diff), width=.2,
                  position=position_dodge(0.05)) +
    theme_bw()+
    theme(axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), panel.border = element_blank(),
          legend.position = "none",
          panel.background = element_blank(), 
          axis.text = element_text(size=18),
          legend.title = element_text(size=12),
          legend.text = element_text(size=12),
          plot.title = element_text(face="bold", size = 20, hjust=0.5),
          plot.subtitle = element_text(face="bold", size = 15, hjust=0.5),
          axis.title=element_text(size=18),
          strip.text = element_text(size=15)) +
    labs(y =  "Difference in the Expected Utility for the Treated",
         x = expression("Cost of Outcome ("*c[0]*")"))
  
  p2 <- ggplot(res, aes(x=c0, y=mean_control_utility_diff)) +
    geom_hline(yintercept = 0,  col = "darkgrey", lty = "dashed") + 
    geom_point()+facet_wrap(~ c1, labeller = label_bquote(atop("Cost of unnecessarily","harsh decision ("*c[1]*")"~"="~.(c1)))) +
    geom_errorbar(aes(ymin=lb_control_utility_diff, ymax=ub_control_utility_diff), width=.2,
                  position=position_dodge(0.05)) +
    theme_bw()+
    theme(axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), panel.border = element_blank(),
          legend.position = "none",
          panel.background = element_blank(), 
          axis.text = element_text(size=18),
          legend.title = element_text(size=12),
          legend.text = element_text(size=12),
          plot.title = element_text(face="bold", size = 20, hjust=0.5),
          plot.subtitle = element_text(face="bold", size = 15, hjust=0.5),
          axis.title=element_text(size=18),
          strip.text = element_text(size=15)) +
    labs(y =  "Difference in the Expected Utility for the Control",
         x = expression("Cost of Outcome ("*c[0]*")"))
  
  return(list(p.treated = p1,
              p.control = p2))
}
