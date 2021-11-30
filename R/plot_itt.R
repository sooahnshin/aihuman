#' Calculate diff-in-means estimates
#' 
#' Calculate average causal effect based on diff-in-means estimator.
#'
#' @param data A \code{data.frame} of which columns includes a binary treatment (Z), an ordinal decision (D), and an outcome variable (Y).
#' 
#' @return A \code{data.frame} of diff-in-means estimates effect for each value of D and Y.
#' 
#' @importFrom stats model.matrix
#' @importFrom stats qnorm
#' 
#' @examples 
#' data(synth)
#' CalDIM(synth)
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
CalDIM <- function(data) {
  sd <- NULL
  
  D = as.factor(data$D)
  D.factor = model.matrix(~0+D)
  dat = as.data.frame(cbind(D.factor, Y = data$Y))
  m = colMeans(dat[data$Z == 1,]) - colMeans(dat[data$Z == 0,])
  s1 = apply(dat[data$Z == 1,], 2, sd)
  s0 = apply(dat[data$Z == 0,], 2, sd)
  
  n0 <- sum(data$Z==0)
  n1 <- sum(data$Z==1)
  s <- sqrt(s0^2/n0+s1^2/n1)
  
  lb <- m + qnorm(0.025)*s
  ub <- m + qnorm(0.975)*s
  
  df = data.frame(mean = m,
                  sd = s,
                  n = n0+n1,
                  lb = lb,
                  ub = ub,
                  itt = colnames(dat))
  
  return(df)
}

#' Calculate diff-in-means estimates
#' 
#' Calculate average causal effect based on diff-in-means estimator.
#'
#' @param data A \code{data.frame} of which columns includes a binary treatment (Z), an ordinal decision (D), and an outcome variable (Y).
#' @param subgroup A list of numeric vectors for the index of each of the five subgroups.
#' @param name.group A character vector including the labels of five subgroups.
#' 
#' @return A \code{data.frame} of diff-in-means estimates for each value of D and Y for each subgroup.
#' 
#' @examples 
#' data(synth)
#' subgroup_synth = list(1:nrow(synth),which(synth$Sex==0),which(synth$Sex==1),
#'                       which(synth$Sex==1&synth$White==0),which(synth$Sex==1&synth$White==1))
#' CalDIMsubgroup(synth, subgroup = subgroup_synth)
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
CalDIMsubgroup <- function(data, 
                   subgroup,
                   name.group = c("Overall", "Female", "Male", "Non-white\nMale", "White\nMale")) {
  # # Age
  # name.group <- c("17~22", "23~28", "29~35", "36~45", "46~")
  
  res  <- CalDIM(data[subgroup[[1]],])
  res <- rbind(res,CalDIM(data[subgroup[[2]],]))
  res <- rbind(res,CalDIM(data[subgroup[[3]],]))
  res <- rbind(res,CalDIM(data[subgroup[[4]],]))
  res <- rbind(res,CalDIM(data[subgroup[[5]],]))
  res$subgroup <- rep(name.group, each = (length(unique(data$D)) + 1))
  
  return(res)
}

#' Plot diff-in-means estimates
#' 
#' See Figure 2 for example.
#'
#' @param res A \code{data.frame} generated with \code{CalDIMsubgroup}.
#' @param y.max Maximum value of y-axis.
#' @param decision.labels Labels of decisions (D).
#' @param col.values Color of point for each decisions.
#' @param shape.values Shape of point for each decisions.
#' 
#' @return A ggplot.
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @import ggplot2
#' 
#' @examples 
#' data(synth)
#' subgroup_synth = list(1:nrow(synth),which(synth$Sex==0),which(synth$Sex==1),
#'                       which(synth$Sex==1&synth$White==0),which(synth$Sex==1&synth$White==1))
#' res_dec = CalDIMsubgroup(synth, subgroup = subgroup_synth)
#' PlotDIMdecisions(res_dec, decision.labels = c("signature","small cash","middle cash","large cash"),
#'                  col.values = c("grey60", "grey30", "grey6", "grey1"), 
#'                  shape.values = c(16, 17, 15, 18))
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
PlotDIMdecisions <- function(res,
                             y.max = 0.2,
                             decision.labels = c("signature bond   ","small cash bond   ","large cash bond"),
                             col.values = c("grey60", "grey30", "grey6"),
                             shape.values = c(16, 17, 15)) {
  itt <- idx <- rgb <- lb <- ub <- NULL
  
  k = nrow(res)/5 - 2
  
  res <- res %>% filter(itt != "Y") 
  res <- res %>%
    mutate(idx = 1:nrow(res)+rep(0:4,each = (k+1)))

  
  p <- ggplot(res, aes(x=idx, y=mean, col = as.factor(itt))) +
    geom_hline(yintercept = 0,  col = "darkgrey", lty = "dashed") +
    geom_vline(xintercept = c(1:max(res$idx))[-res$idx], col = rgb(220,220,220, max = 255)) + 
    geom_point(aes(shape = itt), size = 5) + theme_bw() + labs(x=NULL,y=NULL) +
    geom_errorbar(aes(ymin=lb, ymax=ub), width=.2) + 
    scale_color_manual(name = NULL, values=col.values,
                       labels = decision.labels) +
    scale_shape_manual(name = NULL, values = shape.values,
                       labels = decision.labels) +
    theme(legend.position=c(0.5, 0.09), 
          panel.grid = element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x = element_text(face="bold"),
          axis.text = element_text(size = 20),
          legend.text = element_text(size=20),
          plot.title = element_text(size = 25, face = 'bold', hjust = 0.5),
          legend.direction = "horizontal",
          legend.background = element_rect(fill="transparent",colour=NA)) +
    ylim(-y.max,y.max) + 
    scale_x_continuous(breaks = c(0:4)*(k+2)+(k+2)/2, minor_breaks = NULL,
                       labels = unique(res$subgroup)) + 
    labs(title = "Estimated Average Effects on Judge's Decision")
  
  return(p)
}

#' Plot diff-in-means estimates
#' 
#' See Figure 2 for example.
#'
#' @param res.fta A \code{data.frame} generated with \code{CalDIMsubgroup} with Y = FTA.
#' @param res.nca A \code{data.frame} generated with \code{CalDIMsubgroup} with Y = NCA.
#' @param res.nvca A \code{data.frame} generated with \code{CalDIMsubgroup} with Y = NVCA.
#' @param label.position The position of labels.
#' @param top.margin Top margin of labels.
#' @param bottom.margin Bottom margin of labels.
#' @param y.max Maximum value of y-axis.
#' @param label.size Size of label.
#' @param name.group A character vector including the labels of five subgroups.
#' 
#' @return A ggplot.
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @import ggplot2
#' 
#' @examples 
#' data(synth)
#' subgroup_synth = list(1:nrow(synth),which(synth$Sex==0),which(synth$Sex==1),
#'                       which(synth$Sex==1&synth$White==0),which(synth$Sex==1&synth$White==1))
#' synth_fta <- synth_nca <- synth_nvca <- synth
#' set.seed(123)
#' synth_fta$Y <- sample(0:1, 1000, replace = TRUE)
#' synth_nca$Y <- sample(0:1, 1000, replace = TRUE)
#' synth_nvca$Y <- sample(0:1, 1000, replace = TRUE)
#' res_fta = CalDIMsubgroup(synth_fta, subgroup = subgroup_synth)
#' res_nca = CalDIMsubgroup(synth_nca, subgroup = subgroup_synth)
#' res_nvca = CalDIMsubgroup(synth_nvca, subgroup = subgroup_synth)
#' PlotDIMoutcomes(res_fta, res_nca, res_nvca)
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
PlotDIMoutcomes <- function(res.fta,
                            res.nca,
                            res.nvca,
                            label.position = c("top", "top", "top"), 
                            top.margin = 0.01, bottom.margin = 0.01, 
                            y.max = 0.2,
                            label.size = 7,
                            name.group = c("Overall", "Female", "Male", "Non-white\nMale", "White\nMale")) {
  subgroup <- itt <- idx <- rgb <- lb <- ub <- NULL
  
  res = rbind(res.fta, res.nca, res.nvca)
  res = res[res$itt == "Y",]
  res$itt = rep(c("FTA", "NCA", "NVCA"), each = 5)
  
  res = res %>% arrange(match(subgroup, name.group))
  
  res <- res %>%
    mutate(idx = 1:nrow(res)+rep(0:4,each = 3))
  
  y.position = rep(0, 3)
  for (i in 1:3) {
    if(label.position[i] == "top") {
      y.position[i] = res$ub[i]+top.margin
    } else {
      y.position[i] = res$lb[i]-bottom.margin
    }
  }
  
  p <- ggplot(res, aes(x=idx, y=mean)) +
    geom_hline(yintercept = 0,  col = "darkgrey", lty = "dashed") +
    geom_vline(xintercept = c(1:max(res$idx))[-res$idx], col = rgb(220,220,220, max = 255)) + 
    geom_point(aes(shape = itt), size = 5) + theme_bw() + labs(x=NULL,y=NULL) +
    geom_errorbar(aes(ymin=lb, ymax=ub), width=.2) + 
    scale_shape_manual(name = NULL, values = c(1,2,0),
                       labels = c("FTA   ","NCA   ","NVCA")) +
    theme(legend.position=c(0.5, 0.09), 
          panel.grid = element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x = element_text(face="bold"),
          axis.text = element_text(size = 20),
          legend.text = element_text(size=20),
          plot.title = element_text(size = 25, face = 'bold', hjust = 0.5),
          legend.direction = "horizontal",
          legend.background = element_rect(fill="transparent",colour=NA)) +
    ylim(-y.max,y.max) + 
    scale_x_continuous(breaks = c(0:4)*(2+2)+(2+2)/2, minor_breaks = NULL,
                       labels = unique(res$subgroup)) + 
    annotate(geom="text", x=1, y=y.position[1], label="FTA", size = label.size) + 
    annotate(geom="text", x=2, y=y.position[2], label="NCA", size = label.size) + 
    annotate(geom="text", x=3, y=y.position[3], label="  NVCA", size = label.size) +
    labs(title = "Estimated Average Effect on Arrestee's Behavior")
  
  return(p)
}
