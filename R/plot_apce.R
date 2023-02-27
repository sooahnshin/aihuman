#' Plot APCE
#' 
#' See Figure 4 for example.
#'
#' @param res A \code{data.frame} generated with \code{APCEsummary()}.
#' @param y.max Maximum value of y-axis.
#' @param decision.labels Labels of decisions (D).
#' @param shape.values Shape of point for each decisions.
#' @param col.values Color of point for each principal stratum.
#' @param label A logical argument whether to specify label of each principal stratum. The default is \code{TRUE}.
#' @param r.labels Label of each principal stratum.
#' @param label.position The position of labels.
#' @param top.margin Top margin of labels.
#' @param bottom.margin Bottom margin of labels.
#' @param name.group A character vector including the labels of five subgroups.
#' @param label.size Size of label.
#' 
#' @return A ggplot.
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @import ggplot2
#' 
#' @examples 
#' \donttest{
#' data(synth)
#' sample_mcmc = AiEvalmcmc(data = synth, n.mcmc = 10)
#' subgroup_synth = list(1:nrow(synth),which(synth$Sex==0),which(synth$Sex==1),
#'                       which(synth$Sex==1&synth$White==0),which(synth$Sex==1&synth$White==1))
#' sample_apce = CalAPCE(data = synth, mcmc.re = sample_mcmc, 
#'                       subgroup = subgroup_synth)
#' sample_apce_summary = APCEsummary(sample_apce[["APCE.mcmc"]])
#' PlotAPCE(sample_apce_summary, y.max = 0.25, decision.labels = c("signature","small cash",
#'          "middle cash","large cash"), shape.values = c(16, 17, 15, 18), 
#'          col.values = c("blue", "black", "red", "brown", "purple"), label = FALSE)
#' }
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
PlotAPCE <- function(res, 
                    y.max = 0.1,
                    decision.labels = c("signature bond","small cash bond","large cash bond"),
                    shape.values = c(16, 17, 15),
                    col.values = c("blue", "black", "red", "brown"),
                    label = TRUE,
                    r.labels = c("safe", "easily\npreventable", "prevent-\nable", "risky\n"),
                    label.position = c("top", "top", "top", "top"), 
                    top.margin = 0.01, bottom.margin = 0.01,
                    name.group = c("Overall", "Female", "Male", "Non-white\nMale", "White\nMale"),
                    label.size = 4) {
  subgroup <- Mean <- R <- rgb <- D <- `2.5%` <- `97.5%` <- NULL
  
  for (i in 1:5) {
    res$subgroup[res$subgroup==unique(res$subgroup)[i]] = name.group[i]
  }
  res <- res %>%
    arrange(factor(subgroup, levels = name.group)) 
  
  k = length(unique(res$D)) - 1
  
  res <- res %>%
    mutate(idx = 1:nrow(res) + rep(0:((k+2)*5-1),each = (k+1)) + rep(k*(0:4),each = (k+1)*(k+2))) 
  
  idx <- res %>%
    group_by(subgroup) %>%
    summarise(min = min(idx), max = max(idx)) %>%
    arrange(factor(subgroup, levels = name.group)) 
  
  x.bar = (idx$min[2:5] + idx$max[1:4])/2
  
  res$D = as.factor(res$D)
  res$R = as.factor(res$R)
  

  p <- ggplot(res, aes(x=idx, y=Mean, col = R)) +
    geom_hline(yintercept = 0,  col = "darkgrey", lty = "dashed") +
    geom_vline(xintercept = x.bar, col = rgb(220,220,220, max = 255)) +
    geom_point(aes(shape = D), size = 3) + theme_bw() + labs(x=NULL,y=NULL) +
    geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.2) + 
    scale_shape_manual(name = NULL, values = shape.values,
                       labels = decision.labels) +
    scale_color_manual(values=col.values) +
    theme(legend.position = 'bottom', panel.grid = element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x = element_text(face="bold"),
          axis.text = element_text(size=12),
          legend.text = element_text(size=12),
          plot.title = element_text(face="bold", size = 15, hjust=0.5)) +
    ylim(-y.max,y.max) + 
    scale_x_continuous(breaks = (idx$min + idx$max)/2, 
                       minor_breaks = NULL,
                       labels = name.group) +
    guides(colour="none")
  
  if (label) {
    y.position = rep(0, (k+2))
    ub = res %>%
      filter(subgroup == name.group[1]) %>%
      group_by(R) %>%
      summarise(ub = max(`97.5%`)) %>%
      pull(ub)
    lb = res %>%
      filter(subgroup == name.group[1]) %>%
      group_by(R) %>%
      summarise(lb = min(`2.5%`)) %>%
      pull(lb)
    for (i in 1:(k+2)) {
      if(label.position[i] == "top") {
        y.position[i] = ub[i]+top.margin
      } else {
        y.position[i] = lb[i]-bottom.margin
      }
      p <- p + 
        annotate(geom="text", x=(i-1)*(k+2) + (k+2)/2, y=y.position[i], label=r.labels[i], size = label.size, col = col.values[i])
    }
  }
  
  return(p)
}
