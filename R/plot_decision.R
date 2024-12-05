#' Stacked barplot for the distribution of the decision given psa
#'
#' See Figure 1 for example.
#'
#' @param data A \code{data.frame} of which columns includes an ordinal decision (D), and psa variables (fta, nca, and nvca).
#' @param fta.label Column name of fta score in the data. The default is \code{"FTAScore"}.
#' @param nca.label Column name of nca score in the data. The default is \code{"NCAScore"}.
#' @param nvca.label Column name of nvca score in the data. The default is \code{"NVCAFlag"}.
#' @param d.colors The color of each decision.
#' @param d.labels The label of each decision.
#' @param legend.position The position of legend. The default is \code{"none"}.
#'
#' @return A list of three ggplots.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr summarise
#' @importFrom dplyr n
#' @importFrom dplyr select
#' @importFrom dplyr lag
#' @importFrom dplyr arrange
#' @importFrom dplyr left_join
#' @import ggplot2
#'
#' @examples
#' data(psa_synth)
#' # Control group (PSA not provided)
#' PlotStackedBar(psa_synth[psa_synth$Z == 0, ], d.colors = c(
#'   "grey80", "grey60",
#'   "grey30", "grey10"
#' ), d.labels = c(
#'   "signature", "small",
#'   "middle", "large"
#' ))
#' # Treated group (PSA provided)
#' PlotStackedBar(psa_synth[psa_synth$Z == 0, ], d.colors = c(
#'   "grey80", "grey60",
#'   "grey30", "grey10"
#' ), d.labels = c(
#'   "signature", "small",
#'   "middle", "large"
#' ))
#'
#' @useDynLib aihuman, .registration=TRUE
#' @export
#'
PlotStackedBar <- function(data,
                           fta.label = "FTAScore",
                           nca.label = "NCAScore",
                           nvca.label = "NVCAFlag",
                           d.colors = c("grey60", "grey30", "grey10"),
                           d.labels = c("signature bond", "small cash bond", "large cash bond"),
                           legend.position = "none") {
  FTAScore <- NCAScore <- NVCAFlag <- D <- prop <- yh <- prop2 <- xh <- xl <- yl <- NULL

  colnames(data)[colnames(data) == fta.label] <- "FTAScore"
  colnames(data)[colnames(data) == nca.label] <- "NCAScore"
  colnames(data)[colnames(data) == nvca.label] <- "NVCAFlag"

  dat1 <- as.data.frame(table(data$FTAScore, data$D), stringsAsFactors = F)
  colnames(dat1) <- c("FTAScore", "D", "n")
  dat1$FTAScore <- as.numeric(dat1$FTAScore)
  dat1$D <- as.numeric(dat1$D)
  dat1 <- dat1 %>%
    group_by(FTAScore) %>%
    mutate(prop = 100 * n / sum(n)) %>%
    left_join(data %>%
      group_by(FTAScore) %>%
      summarise(n = n()) %>%
      mutate(prop2 = 100 * n / sum(n)) %>%
      select(-n)) %>%
    arrange(FTAScore, -D) %>%
    mutate(yh = cumsum(prop)) %>%
    mutate(yl = lag(yh, default = 0)) %>%
    group_by(D) %>%
    mutate(xh = cumsum(prop2)) %>%
    mutate(xl = lag(xh, default = 0))

  p1 <- dat1 %>%
    ggplot(aes(
      xmin = xl + FTAScore * 2, xmax = xh + FTAScore * 2,
      ymin = yl, ymax = yh,
      group = as.factor(D), fill = as.factor(D)
    )) +
    geom_rect() +
    scale_fill_manual(
      values = d.colors,
      labels = d.labels,
      name = "Judge's Decision"
    ) +
    theme_bw() +
    theme(
      legend.position = legend.position,
      panel.grid = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(face = "bold"),
      axis.text = element_text(size = 20),
      plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20)
    ) +
    scale_x_continuous(
      breaks = round(unique((dat1$xh - dat1$xl) / 2 + dat1$xl + dat1$FTAScore * 2)),
      labels = unique(dat1$FTAScore)
    ) +
    scale_y_continuous(
      breaks = c(0, 25, 50, 75, 100),
      labels = paste0(c(0, 25, 50, 75, 100), "%")
    ) +
    labs(title = "Failure to Appear\n(FTA)")

  dat2 <- as.data.frame(table(data$NCAScore, data$D), stringsAsFactors = F)
  colnames(dat2) <- c("NCAScore", "D", "n")
  dat2$NCAScore <- as.numeric(dat2$NCAScore)
  dat2$D <- as.numeric(dat2$D)
  dat2 <- dat2 %>%
    group_by(NCAScore) %>%
    mutate(prop = 100 * n / sum(n)) %>%
    left_join(data %>%
      group_by(NCAScore) %>%
      summarise(n = n()) %>%
      mutate(prop2 = 100 * n / sum(n)) %>%
      select(-n)) %>%
    arrange(NCAScore, -D) %>%
    mutate(yh = cumsum(prop)) %>%
    mutate(yl = lag(yh, default = 0)) %>%
    group_by(D) %>%
    mutate(xh = cumsum(prop2)) %>%
    mutate(xl = lag(xh, default = 0))

  p2 <- dat2 %>%
    ggplot(aes(
      xmin = xl + NCAScore * 2, xmax = xh + NCAScore * 2,
      ymin = yl, ymax = yh,
      group = as.factor(D), fill = as.factor(D)
    )) +
    geom_rect() +
    scale_fill_manual(
      values = d.colors,
      labels = d.labels,
      name = "Judge's Decision"
    ) +
    theme_bw() +
    theme(
      legend.position = legend.position,
      panel.grid = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(face = "bold"),
      axis.text = element_text(size = 20),
      plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20)
    ) +
    scale_x_continuous(
      breaks = round(unique((dat2$xh - dat2$xl) / 2 + dat2$xl + dat2$NCAScore * 2)),
      labels = unique(dat2$NCAScore)
    ) +
    scale_y_continuous(
      breaks = c(0, 25, 50, 75, 100),
      labels = paste0(c(0, 25, 50, 75, 100), "%")
    ) +
    labs(title = "New Criminal Activity\n(NCA)")

  dat3 <- as.data.frame(table(data$NVCAFlag, data$D), stringsAsFactors = F)
  colnames(dat3) <- c("NVCAFlag", "D", "n")
  dat3$NVCAFlag <- as.numeric(dat3$NVCAFlag)
  dat3$D <- as.numeric(dat3$D)
  dat3 <- dat3 %>%
    group_by(NVCAFlag) %>%
    mutate(prop = 100 * n / sum(n)) %>%
    left_join(data %>%
      group_by(NVCAFlag) %>%
      summarise(n = n()) %>%
      mutate(prop2 = 100 * n / sum(n)) %>%
      select(-n)) %>%
    arrange(NVCAFlag, -D) %>%
    mutate(yh = cumsum(prop)) %>%
    mutate(yl = lag(yh, default = 0)) %>%
    group_by(D) %>%
    mutate(xh = cumsum(prop2)) %>%
    mutate(xl = lag(xh, default = 0))

  p3 <- dat3 %>%
    ggplot(aes(
      xmin = xl + NVCAFlag * 2, xmax = xh + NVCAFlag * 2,
      ymin = yl, ymax = yh,
      group = as.factor(D), fill = as.factor(D)
    )) +
    geom_rect() +
    scale_fill_manual(
      values = d.colors,
      labels = d.labels,
      name = "Judge's Decision"
    ) +
    theme_bw() +
    theme(
      legend.position = legend.position,
      panel.grid = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(face = "bold"),
      axis.text = element_text(size = 20),
      plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      legend.text.align = 0
    ) +
    scale_x_continuous(
      breaks = round(unique((dat3$xh - dat3$xl) / 2 + dat3$xl + dat3$NVCAFlag * 2)),
      labels = unique(dat3$NVCAFlag)
    ) +
    scale_y_continuous(
      breaks = c(0, 25, 50, 75, 100),
      labels = paste0(c(0, 25, 50, 75, 100), "%")
    ) +
    labs(title = "New Violent Criminal Activity\n(NVCA)")
  return(list(p.fta = p1, p.nca = p2, p.nvca = p3))
}

#' Stacked barplot for the distribution of the decision given DMF recommendation
#'
#' See Figure 1 for example.
#'
#' @param data A \code{data.frame} of which columns includes a binary treatment (Z; PSA provision), an ordinal decision (D), and DMF recommendation.
#' @param dmf.label Column name of DMF recommendation in the data. The default is \code{"dmf"}.
#' @param dmf.category The name of each category of DMF recommendation.
#' @param d.colors The color of each decision.
#' @param d.labels The label of each decision.
#' @param legend.position The position of legend. The default is \code{"none"}.
#'
#' @return A list of three ggplots.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr summarise
#' @importFrom dplyr n
#' @importFrom dplyr select
#' @importFrom dplyr lag
#' @importFrom dplyr arrange
#' @importFrom dplyr left_join
#' @import ggplot2
#'
#' @examples
#' data(psa_synth)
#' PlotStackedBarDMF(psa_synth, dmf.label = "DMF", d.colors = c(
#'   "grey80",
#'   "grey60", "grey30", "grey10"
#' ), d.labels = c(
#'   "signature",
#'   "small", "middle", "large"
#' ))
#'
#' @useDynLib aihuman, .registration=TRUE
#' @export
#'
PlotStackedBarDMF <- function(data,
                              dmf.label = "dmf",
                              dmf.category = NULL,
                              d.colors = c("grey60", "grey30", "grey10"),
                              d.labels = c("signature bond", "small cash bond", "large cash bond"),
                              legend.position = "none") {
  dmf <- D <- prop <- yh <- prop2 <- xh <- xl <- yl <- NULL

  colnames(data)[colnames(data) == dmf.label] <- "dmf"

  ctrl <- data[data$Z == 0, ]
  treat <- data[data$Z == 1, ]

  dat1 <- as.data.frame(table(treat$dmf, treat$D), stringsAsFactors = F)
  colnames(dat1) <- c("dmf", "D", "n")
  dat1$dmf <- as.numeric(dat1$dmf)
  dat1$D <- as.numeric(dat1$D)
  dat1 <- dat1 %>%
    group_by(dmf) %>%
    mutate(prop = 100 * n / sum(n)) %>%
    left_join(treat %>%
      group_by(dmf) %>%
      summarise(n = n()) %>%
      mutate(prop2 = 100 * n / sum(n)) %>%
      select(-n)) %>%
    arrange(dmf, -D) %>%
    mutate(yh = cumsum(prop)) %>%
    mutate(yl = lag(yh, default = 0)) %>%
    group_by(D) %>%
    mutate(xh = cumsum(prop2)) %>%
    mutate(xl = lag(xh, default = 0))
  if (is.null(dmf.category)) dmf.category <- unique(dat1$dmf)
  p1 <- dat1 %>%
    ggplot(aes(
      xmin = xl + dmf * 2, xmax = xh + dmf * 2,
      ymin = yl, ymax = yh,
      group = as.factor(D), fill = as.factor(D)
    )) +
    geom_rect() +
    scale_fill_manual(
      values = d.colors,
      labels = d.labels,
      name = "Judge's Decision"
    ) +
    theme_bw() +
    theme(
      legend.position = legend.position,
      panel.grid = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(face = "bold"),
      axis.text = element_text(size = 20),
      plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      legend.text.align = 0
    ) +
    scale_x_continuous(
      breaks = round(unique((dat1$xh - dat1$xl) / 2 + dat1$xl + dat1$dmf * 2)),
      labels = dmf.category
    ) +
    scale_y_continuous(
      breaks = c(0, 25, 50, 75, 100),
      labels = paste0(c(0, 25, 50, 75, 100), "%")
    ) +
    labs(title = "DMF\nRecommendation")

  dat2 <- as.data.frame(table(ctrl$dmf, ctrl$D), stringsAsFactors = F)
  colnames(dat2) <- c("dmf", "D", "n")
  dat2$dmf <- as.numeric(dat2$dmf)
  dat2$D <- as.numeric(dat2$D)
  dat2 <- dat2 %>%
    group_by(dmf) %>%
    mutate(prop = 100 * n / sum(n)) %>%
    left_join(ctrl %>%
      group_by(dmf) %>%
      summarise(n = n()) %>%
      mutate(prop2 = 100 * n / sum(n)) %>%
      select(-n)) %>%
    arrange(dmf, -D) %>%
    mutate(yh = cumsum(prop)) %>%
    mutate(yl = lag(yh, default = 0)) %>%
    group_by(D) %>%
    mutate(xh = cumsum(prop2)) %>%
    mutate(xl = lag(xh, default = 0))
  if (is.null(dmf.category)) dmf.category <- unique(dat2$dmf)
  p2 <- dat2 %>%
    ggplot(aes(
      xmin = xl + dmf * 2, xmax = xh + dmf * 2,
      ymin = yl, ymax = yh,
      group = as.factor(D), fill = as.factor(D)
    )) +
    geom_rect() +
    scale_fill_manual(
      values = d.colors,
      labels = d.labels,
      name = "Judge's Decision"
    ) +
    theme_bw() +
    theme(
      legend.position = legend.position,
      panel.grid = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(face = "bold"),
      axis.text = element_text(size = 20),
      plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      legend.text.align = 0
    ) +
    scale_x_continuous(
      breaks = round(unique((dat2$xh - dat2$xl) / 2 + dat2$xl + dat2$dmf * 2)),
      labels = dmf.category
    ) +
    scale_y_continuous(
      breaks = c(0, 25, 50, 75, 100),
      labels = paste0(c(0, 25, 50, 75, 100), "%")
    ) +
    labs(title = "DMF\nRecommendation")

  return(list(p.treat = p1, p.control = p2))
}

#' Pulling ggplot legend
#'
#' @param p A \code{ggplot} of which legend should be pulled.
#'
#' @return A ggplot legend.
#'
#' @import ggplot2
#'
#' @useDynLib aihuman, .registration=TRUE
#' @export
#'
g_legend <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]

  return(legend)
}
