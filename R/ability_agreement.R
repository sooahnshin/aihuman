#' Agreement of Human and AI Decision Makers
#'
#' Estimate the impact of AI recommendations on the agreement between human decisions and AI recommendations using a difference-in-means estimator of an indicator \eqn{1\{D_i = A_i\}}.
#'
#' @param Y An observed outcome (binary: numeric vector of 0 or 1).
#' @param D An observed decision (binary: numeric vector of 0 or 1).
#' @param Z A treatment indicator (binary: numeric vector of 0 or 1).
#' @param A An AI recommendation (binary: numeric vector of 0 or 1).
#' @param X Pretreatment covariate used for subgroup analysis (vector). Must be the same length as Y, D, Z, and A if provided. Default is NULL.
#'
#' @return A tibble with the following columns:
#' \itemize{
#'   \item \code{X}: Pretreatment covariate (if provided).
#'   \item \code{agree_diff}: Difference in agreement between human decisions and AI recommendations.
#'   \item \code{agree_diff_se}: Standard error of the difference in agreement.
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise cross_join filter mutate rename select
#' @examples
#' compute_stats_agreement(
#'   Y = NCAdata$Y,
#'   D = ifelse(NCAdata$D == 0, 0, 1),
#'   Z = NCAdata$Z,
#'   A = PSAdata$DMF
#' )
#'
#' @useDynLib aihuman, .registration = TRUE
#' @export
compute_stats_agreement <- function(Y, D, Z, A, X = NULL) {
  # Check that Y, D, Z, and A are binary
  if (!all(Y %in% c(0, 1))) stop("Y must be binary")
  if (!all(D %in% c(0, 1))) stop("D must be binary")
  if (!all(Z %in% c(0, 1))) stop("Z must be binary")
  if (!all(A %in% c(0, 1))) stop("A must be binary")

  # Check lengths using all.equal
  if (is.null(X)) {
    if (!isTRUE(all.equal(length(Y), length(D), length(Z), length(A)))) {
      stop("Y, D, Z, and A must have the same length")
    }
  } else if (is.vector(X)) {
    if (!isTRUE(all.equal(length(X), length(Y), length(D), length(Z), length(A)))) {
      stop("X must have the same length as Y, D, Z, and A")
    }
  } else {
    stop("X must be either NULL or a vector")
  }

  # Use a placeholder if X is NULL
  flag <- is.null(X)
  if (flag) {
    X <- rep(1, length(Y))
  }

  # Create a data frame and calculate statistics
  props <- data.frame(Y = Y, D = D, Z = Z, A = A, X = X) %>%
    mutate(W = ifelse(D == A, 1, 0)) %>%
    group_by(Z, X) %>%
    summarise(
      Wbar = mean(W),
      sigma_W = sd(W),
      n = n(),
      .groups = "drop"
    )

  # Compute differences in agreement and standard errors
  out <- cross_join(props, props) %>%
    filter(X.y == X.x) %>%
    mutate(
      agree_diff = Wbar.x - Wbar.y,
      agree_diff_se = sqrt(sigma_W.x^2 / n.x + sigma_W.y^2 / n.y)
    ) %>%
    filter(Z.x > Z.y) %>%
    select(X.x, contains("diff")) %>%
    rename(X = X.x)

  # Return the result
  if (flag) {
    return(out %>% select(-X))
  } else {
    return(out)
  }
}

#' Visualize Agreement
#'
#' Visualize the agreement between human decisions and AI recommendations using a difference-in-means estimator of an indicator \eqn{1\{D_i = A_i\}}.
#' Generate a plot based on the overall agreement and subgroup-specific agreement.
#'
#' @param Y An observed outcome (binary: numeric vector of 0 or 1).
#' @param D An observed decision (binary: numeric vector of 0 or 1).
#' @param Z A treatment indicator (binary: numeric vector of 0 or 1).
#' @param A An AI recommendation (binary: numeric vector of 0 or 1).
#' @param subgroup1 A pretreatment covariate used for subgroup analysis (vector).
#' @param subgroup2 A pretreatment covariate used for subgroup analysis (vector).
#' @param label.subgroup1 A label for subgroup1 (character). Default "Subgroup 1".
#' @param label.subgroup2 A label for subgroup2 (character). Default "Subgroup 2".
#' @param x.order An order for the x-axis (character vector). Default NULL.
#' @param p.title A title for the plot (character). Default NULL.
#' @param p.lb A lower bound for the y-axis (numeric). Default -0.3.
#' @param p.ub An upper bound for the y-axis (numeric). Default 0.3.
#' @param y.lab A label for the y-axis (character). Default "Impact of PSA".
#'
#' @return A ggplot object.
#'
#' @importFrom dplyr case_when
#' @importFrom stringr str_detect
#'
#' @examples
#' plot_agreement(
#'   Y = NCAdata$Y,
#'   D = ifelse(NCAdata$D == 0, 0, 1),
#'   Z = NCAdata$Z,
#'   A = PSAdata$DMF,
#'   subgroup1 = ifelse(NCAdata$White == 1, "White", "Non-white"),
#'   subgroup2 = ifelse(NCAdata$Sex == 1, "Male", "Female"),
#'   label.subgroup1 = "Race",
#'   label.subgroup2 = "Gender",
#'   x.order = c("Overall", "Non-white", "White", "Female", "Male")
#' )
#'
#' @useDynLib aihuman, .registration = TRUE
#' @export
plot_agreement <- function(Y, D, Z, A,
                           subgroup1, subgroup2,
                           label.subgroup1 = "Subgroup 1", label.subgroup2 = "Subgroup 2",
                           x.order = NULL,
                           p.title = NULL, p.lb = -0.3, p.ub = 0.3, y.lab = "Impact of PSA") {
  diff_agree <- compute_stats_agreement(Y = Y, D = D, Z = Z, A = A, X = NULL) |>
    mutate(cov = "Overall", X = "Overall") |>
    relocate(cov, X) |>
    bind_rows(compute_stats_agreement(Y = Y, D = D, Z = Z, A = A, X = subgroup1) |> mutate(cov = label.subgroup1) |> relocate(cov)) |>
    bind_rows(compute_stats_agreement(Y = Y, D = D, Z = Z, A = A, X = subgroup2) |> mutate(cov = label.subgroup2) |> relocate(cov))
  diff_vis <- diff_agree |>
    pivot_longer(cols = -c("cov", "X")) |>
    mutate(
      metric = case_when(
        str_detect(name, "agree_diff") ~ "Agreement of Human and PSA Recommendations",
        TRUE ~ ""
      ),
      se = ifelse(str_detect(name, "_se"), "se", "est")
    ) |>
    select(-name) |>
    pivot_wider(values_from = "value", names_from = "se") |>
    mutate(cov = factor(cov, levels = c("Overall", label.subgroup1, label.subgroup2)))
  if (!is.null(x.order)) {
    diff_vis <- diff_vis |>
      mutate(X = factor(X, levels = x.order))
  }
  p_vis <- vis_agreement(df = diff_vis, p.title = p.title, p.lb = p.lb, p.ub = p.ub, y.lab = y.lab)
  return(p_vis)
}

#' Visualize Agreement (internal)
#'
#' Internal function to visualize the agreement between human decisions and AI recommendations using a difference-in-means estimator of an indicator \eqn{1\{D_i = A_i\}}.
#'
#' @param df A data frame generated by \code{compute_stats_agreement} and reshaped.
#' @param p.title A title for the plot (character). Default NULL.
#' @param p.lb A lower bound for the y-axis (numeric). Default -0.2.
#' @param p.ub An upper bound for the y-axis (numeric). Default 0.2.
#' @param y.lab A label for the y-axis (character). Default "Impact of PSA".
#'
#' @return A ggplot object.
#'
#' @keywords internal
#' @useDynLib aihuman, .registration = TRUE
#'
vis_agreement <- function(df, p.title = NULL, p.lb = -0.2, p.ub = 0.2, y.lab = "Impact of PSA") {
  annotations <- data.frame(
    x = c(1 - 0.2, 1 - 0.2),
    y.start = c(p.ub - 0.3 * p.ub, p.lb + 0.3 * p.ub),
    y.end = c(p.ub, p.lb),
    label = c("More agreement", "Less agreement")
  )
  annotations_bottom <- annotations %>%
    filter(label == "More agreement")
  annotations_top <- annotations %>%
    filter(label == "Less agreement")
  df %>%
    ggplot(aes(x = X, color = cov, shape = cov)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_pointrange(aes(y = est, ymin = est - qnorm(0.975) * se, ymax = est + qnorm(0.975) * se),
      position = position_dodge(.5), lineend = "round", linewidth = 1.2, size = 1.2
    ) +
    ylab(y.lab) +
    scale_color_brewer("", type = "qual", palette = "Set1") +
    ggtitle(p.title) +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    scale_y_continuous(limits = c(p.lb, p.ub)) +
    geom_text(
      data = annotations_bottom, aes(x = x, y = y.end, label = label),
      inherit.aes = FALSE, vjust = -1, hjust = 1.2, angle = 90
    ) +
    geom_text(
      data = annotations_top, aes(x = x, y = y.end, label = label),
      inherit.aes = FALSE, vjust = -1, hjust = -0.2, angle = 90
    ) +
    geom_segment(
      data = annotations, aes(x = x, xend = x, y = y.start, yend = y.end),
      arrow = arrow(length = unit(0.1, "inches")), inherit.aes = FALSE
    )
}

#' Table of Agreement
#'
#' Estimate the impact of AI recommendations on the agreement between human decisions and AI recommendations using a difference-in-means estimator of an indicator \eqn{1\{D_i = A_i\}}.
#' Generate a table based on the overall agreement and subgroup-specific agreement.
#'
#' @param Y An observed outcome (binary: numeric vector of 0 or 1).
#' @param D An observed decision (binary: numeric vector of 0 or 1).
#' @param Z A treatment indicator (binary: numeric vector of 0 or 1).
#' @param A An AI recommendation (binary: numeric vector of 0 or 1).
#' @param subgroup1 A pretreatment covariate used for subgroup analysis (vector).
#' @param subgroup2 A pretreatment covariate used for subgroup analysis (vector).
#' @param label.subgroup1 A label for subgroup1 (character). Default "Subgroup 1".
#' @param label.subgroup2 A label for subgroup2 (character). Default "Subgroup 2".
#'
#' @return A tibble with the following columns:
#' \itemize{
#'   \item \code{cov}: Subgroup label.
#'   \item \code{X}: Subgroup value.
#'   \item \code{agree_diff}: Difference in agreement between human decisions and AI recommendations.
#'   \item \code{agree_diff_se}: Standard error of the difference in agreement.
#'   \item \code{ci_lb}: Lower bound of the 95\% confidence interval.
#'   \item \code{ci_ub}: Upper bound of the 95\% confidence interval.
#' }
#' @importFrom dplyr mutate relocate bind_rows
#' @examples
#' table_agreement(
#'   Y = NCAdata$Y,
#'   D = ifelse(NCAdata$D == 0, 0, 1),
#'   Z = NCAdata$Z,
#'   A = PSAdata$DMF,
#'   subgroup1 = ifelse(NCAdata$White == 1, "White", "Non-white"),
#'   subgroup2 = ifelse(NCAdata$Sex == 1, "Male", "Female"),
#'   label.subgroup1 = "Race",
#'   label.subgroup2 = "Gender"
#' )
#'
#' @useDynLib aihuman, .registration = TRUE
#' @export
table_agreement <- function(Y, D, Z, A,
                            subgroup1, subgroup2,
                            label.subgroup1 = "Subgroup 1", label.subgroup2 = "Subgroup 2") {
  diff_agree <- compute_stats_agreement(Y = Y, D = D, Z = Z, A = A, X = NULL) |>
    mutate(cov = "Overall", X = "Overall") |>
    relocate(cov, X) |>
    bind_rows(compute_stats_agreement(Y = Y, D = D, Z = Z, A = A, X = subgroup1) |> mutate(cov = label.subgroup1) |> relocate(cov)) |>
    bind_rows(compute_stats_agreement(Y = Y, D = D, Z = Z, A = A, X = subgroup2) |> mutate(cov = label.subgroup2) |> relocate(cov))
  tab <- diff_agree |>
    mutate(ci_ub = agree_diff + qnorm(0.975) * agree_diff_se, ci_lb = agree_diff - qnorm(0.975) * agree_diff_se)
  return(tab)
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "W", "Wbar.x", "Wbar.y", "X", "X.x", "X.y", "Z.x", "Z.y",
    "agree_diff", "agree_diff_se", "contains", "cov", "est",
    "label", "n.x", "n.y", "name", "sd", "se",
    "sigma_W.x", "sigma_W.y", "x", "y.end", "y.start"
  ))
}
