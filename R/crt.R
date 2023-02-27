#' Conduct conditional randomization test
#' 
#' See S3.1 for more details.
#'
#' @param D A numeric vector of judge's decision.
#' @param Z A numeric vector of treatment variable.
#' @param CourtEvent_HearingDate The court event hearing date.
#' @param n Number of permutations.
#' @param seed.number An integer for random number generator.
#' 
#' @return A \code{list} of the observed and permuted test statistics and its p-value.
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate
#' @importFrom dpltr left_join
#' @importFrom dpltr select
#' @importFrom dpltr filter
#' @importFrom purrr rbernoulli
#' @importFrom MASS polr
#' 
#' @examples 
#' \dontrun{
#' data(synth)
#' data(hearingdate_synth)
#' crt <- SpilloverCRT(D = synth$D, Z = synth$Z, CourtEvent_HearingDate = hearingdate_synth)
#' }
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
SpilloverCRT <- function(D, Z, CourtEvent_HearingDate,
                         n = 100, seed.number = 12345) {
  lm <- rand.Z <- .Ztilde <- NULL
  
  dat = data.frame(D = D,
                   Z = Z,
                   CourtEvent_HearingDate = as.factor(CourtEvent_HearingDate))
  dat.Ztilde <- dat %>%
    group_by(CourtEvent_HearingDate) %>%
    summarise(.Ztilde = mean(Z==1), .groups = 'drop') %>% # .Ztilde = proportion of treated per date
    mutate(Ztilde = lag(.Ztilde)) # Ztilde = proportion of treated for one previous date
  
  dat <- dat %>%
    left_join(dat.Ztilde, by = "CourtEvent_HearingDate") 
  
  n_dat = nrow(dat)
  even.days = levels(dat$CourtEvent_HearingDate)[seq(2,length(levels(dat$CourtEvent_HearingDate)),2)]
  
  dat.even <- dat %>% filter(CourtEvent_HearingDate %in% even.days)
  t.stat = rep(NA, n) # test statistic
  mod <- lm(D ~ Z + Ztilde, data = dat.even)
  t.obs <- (mod$coefficients["Ztilde"])^2 # observed statstic = omegaZtilde.k
  
  set.seed(seed.number)
  for (s in 1:n) {
    dat.Ztilde.s <- dat %>%
      mutate(rand.Z = rbernoulli(n(), 0.5)) %>%
      group_by(CourtEvent_HearingDate) %>%
      summarise(.Ztilde = mean(rand.Z==1), .groups = 'drop') %>% # .Ztilde = proportion of treated per date
      mutate(Ztilde = lag(.Ztilde)) # Ztilde = proportion of treated for one previous date
    dat.s <- dat %>%
      select(D, Z, CourtEvent_HearingDate) %>%
      left_join(dat.Ztilde.s, by = "CourtEvent_HearingDate") %>%
      filter(CourtEvent_HearingDate %in% even.days)
    mod <- polr(D ~ Z + Ztilde, data = dat.s %>% mutate(D = as.factor(D)))
    t.stat[s] = (mod$coefficients["Ztilde"])^2
    # cat(s,"th done!\n")
  }
  save_pval = mean(t.stat>=t.obs)
  
  return(list(permuted.statistics = t.stat,
              observed.statistics = t.obs,
              pval = save_pval))
}

#' Plot conditional randomization test
#' 
#' See Figure S8 for example.
#'
#' @param res A \code{list} generated with \code{SpilloverCRT}.
#' 
#' @return A ggplot
#' 
#' @import ggplot2
#' 
#' @examples 
#' \dontrun{
#' data(synth)
#' data(hearingdate_synth)
#' crt <- SpilloverCRT(D = synth$D, Z = synth$Z, CourtEvent_HearingDate = hearingdate_synth)
#' PlotSpilloverCRT(crt)
#' }
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
PlotSpilloverCRT <- function(res) {
  `..density..` <- density <- NULL
  
  t.stat = res$permuted.statistics
  t.obs = res$observed.statistics
  save_pval = res$pval
  
  p <- ggplot(data.frame(t = t.stat), aes(x=t)) + 
    geom_histogram(binwidth=0.005,
                   color="black", fill="white",
                   aes(y=..density..)) +
    geom_density() +
    theme_bw() +
    xlab(expression(T)) +
    ylab("Density")+
    theme(legend.position="none", 
          panel.grid = element_blank(),
          axis.ticks.x=element_blank(),
          axis.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          axis.title = element_text(size = 20)) +
    geom_vline(xintercept = t.obs, color = "red") +
    annotate("text", x = t.obs, y = max(density(t.stat)$y), 
             label = paste("p-value = ", round(save_pval, 4)),
             hjust=-0.1, size = 6)
  
  return(p)
}
#' Conduct power analysis of conditional randomization test
#' 
#' See S3.2 for more details.
#'
#' @param D A numeric vector of judge's decision.
#' @param Z A numeric vector of treatment variable.
#' @param CourtEvent_HearingDate The court event hearing date.
#' @param n Number of permutations.
#' @param m Number of permutations.
#' @param size The number of parallel computing. The default is \code{2}.
#' @param cand_omegaZtilde Candidate values
#' 
#' @return A \code{data.frame} of the result of power analysis.
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate
#' @importFrom dpltr left_join
#' @importFrom dpltr select
#' @importFrom dpltr filter
#' @importFrom purrr rbernoulli
#' @importFrom foreach foreach
#' @importFrom foreach %do%
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom MASS polr
#' 
#' @examples 
#' \dontrun{
#' data(synth)
#' data(hearingdate_synth)
#' crt_power <- SpilloverCRTpower(D = synth$D, Z = synth$Z, 
#'                                CourtEvent_HearingDate = hearingdate_synth,
#'                                size = 1) # adjust the size
#' }
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
SpilloverCRTpower <- function(D, Z, CourtEvent_HearingDate,
                              n = 4, m = 4, size = 2,
                              cand_omegaZtilde = seq(-1.5, 1.5, by = 0.5)) {
  if(size == 1) {
    message("Increase the size.")
  } else {
    omega <- prop <- predict <- m_lb <- m_ub <- lm <- rand.Z <- .Ztilde <- NULL
    
    numCores <- detectCores()
    registerDoParallel(numCores)
    size.cores = floor(numCores/size)
    
    dat = data.frame(D = D,
                     Z = Z,
                     CourtEvent_HearingDate = as.factor(CourtEvent_HearingDate))
    
    dat.Ztilde <- dat %>%
      group_by(CourtEvent_HearingDate) %>%
      summarise(.Ztilde = mean(Z==1), .groups = 'drop') %>% # .Ztilde = proportion of treated per date
      mutate(Ztilde = lag(.Ztilde)) # Ztilde = proportion of treated for one previous date
    
    dat <- dat %>%
      left_join(dat.Ztilde, by = "CourtEvent_HearingDate") 
    
    n_dat = nrow(dat)
    even.days = levels(dat$CourtEvent_HearingDate)[seq(2,length(levels(dat$CourtEvent_HearingDate)),2)]
    
    ### Step 1-1: Fit the model to the observed data
    mod.ord.logit <- polr(D ~ Z + Ztilde, data = dat %>% mutate(D = as.factor(D)))
    
    ### Step 1-2. Set the candidate omega_{Z_tilde}
    prop_ls <- omega_ls <- rep(NA, length(cand_omegaZtilde))
    for (k in 1:length(cand_omegaZtilde)) {
      lb = c(1, floor(m/size.cores)*1:(size.cores-1) + 1)
      ub = c(floor(m/size.cores)*1:(size.cores-1), m)
      
      omegaZtilde.k = cand_omegaZtilde[k]
      mod.ord.logit$coefficients["Ztilde"] = omegaZtilde.k
      Phat.k <- predict(mod.ord.logit, newdata = dat, type="probs")
      save_pval = rep(NA, m)
      foreach (m_lb=lb, m_ub=ub) %do% {
        for (j in m_lb:m_ub) {
          set.seed(j)
          D.k <- rep(NA, n_dat)
          ### Step 2: Generate D according to the model
          for (i in (1:n_dat)[!is.na(dat$Ztilde)]) {
            D.k[i] <- sample(
              x = unique(D), 
              size = 1, 
              prob = Phat.k[i,]
            )
          }
          
          dat.k <- data.frame(D = D.k,
                              Z = dat$Z,
                              Ztilde = dat$Ztilde,
                              CourtEvent_HearingDate = dat$CourtEvent_HearingDate)
          
          ### Step 3: Conduct the test and save the p-value
          dat.k.even <- dat.k %>% filter(CourtEvent_HearingDate %in% even.days)
          # summary(dat.k.even)
          t.stat = rep(NA, n) # test statistic
          mod <- lm(D ~ Z + Ztilde, data = dat.k.even)
          t.obs <- (mod$coefficients["Ztilde"])^2 # observed statstic = omegaZtilde.k
          
          for (s in 1:n) {
            dat.k.Ztilde.s <- dat.k %>%
              mutate(rand.Z = rbernoulli(n(), 0.5)) %>%
              group_by(CourtEvent_HearingDate) %>%
              summarise(.Ztilde = mean(rand.Z==1), .groups = 'drop') %>% # .Ztilde = proportion of treated per date
              mutate(Ztilde = lag(.Ztilde)) # Ztilde = proportion of treated for one previous date
            dat.k.s <- dat.k %>%
              select(D, Z, CourtEvent_HearingDate) %>%
              left_join(dat.k.Ztilde.s, by = "CourtEvent_HearingDate") %>%
              filter(CourtEvent_HearingDate %in% even.days)
            mod <- lm(D ~ Z + Ztilde, data = dat.k.s)
            t.stat[s] = (mod$coefficients["Ztilde"])^2
          }
          save_pval[j] = mean(t.stat>=t.obs)
          # cat(j,"th done!\n")
        }
      }
      prop_ls[k] <- mean(save_pval < 0.05)
      omega_ls[k] <- omegaZtilde.k
    }
    
    pow_df <- data.frame(prop = prop_ls,
                         omega = omega_ls) 
    return(pow_df)
  }
}
#' Plot power analysis of conditional randomization test
#' 
#' See Figure S9 for example.
#'
#' @param res A \code{data.frame} generated with \code{SpilloverCRTpower}.
#' 
#' @return A ggplot
#' 
#' @import ggplot2
#' 
#' @examples 
#' \dontrun{
#' data(synth)
#' data(hearingdate_synth)
#' crt_power <- SpilloverCRTpower(D = synth$D, Z = synth$Z, 
#'                                CourtEvent_HearingDate = hearingdate_synth,
#'                                size = 1) # adjust the size
#' # PlotSpilloverCRTpower(crt_power)
#' }
#' 
#' @useDynLib aihuman, .registration=TRUE
#' @export
#' 
PlotSpilloverCRTpower <- function(res) {
  omega <- prop <- NULL
  
  p <- ggplot(res) +
    geom_point(aes(x = omega, y = prop)) +
    geom_line(aes(x = omega, y = prop)) +
    theme_bw() +
    theme(legend.position="none", 
          axis.ticks.x=element_blank(),
          axis.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          axis.title = element_text(size = 20)) +
    scale_x_continuous(breaks = res$omega) +
    scale_y_continuous(limits = c(0,1)) +
    labs(x = expression(omega),
         y = "Proportion")
  return(p)
}
