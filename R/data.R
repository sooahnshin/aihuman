#' Synthetic data
#'
#' A synthetic dataset containing pre-treatment covariates, a binary treatment (Z), an ordinal decision (D),
#' and an outcome variable (Y).
#'
#' @format A data frame with 1000 rows and 11 variables:
#' \describe{
#'   \item{Z}{binary treatment}
#'   \item{D}{ordinal decision}
#'   \item{Y}{outcome}
#'   \item{Sex}{male or female}
#'   \item{White}{white or non-white}
#'   \item{Age}{age}
#'   \item{CurrentViolentOffense}{binary variable for current violent offense}
#'   \item{PendingChargeAtTimeOfOffense}{binary variable for pending charge (felony, misdemeanor, or both) at the time of offense}
#'   \item{PriorMisdemeanorConviction}{binary variable for prior conviction of misdemeanor}
#'   \item{PriorFelonyConviction}{binary variable for prior conviction of felony}
#'   \item{PriorViolentConviction}{four-level ordinal variable for prior violent conviction}
#' }
"synth"

#' Synthetic PSA data
#'
#' A synthetic dataset containing a binary treatment (Z), ordinal decision (D), three PSA variables (FTAScore, NCAScore, and NVCAFlag), and DMF recommendation.
#'
#' @format A data frame with 1000 rows and 4 variables:
#' \describe{
#'   \item{Z}{binary treatment}
#'   \item{D}{ordinal decision}
#'   \item{FTAScore}{FTA score}
#'   \item{NCAScore}{NCA score}
#'   \item{NVCAFlag}{NVCA flag}
#'   \item{DMF}{DMF recommendation}
#' }
"psa_synth"

#' Synthetic court event hearing date
#'
#' A synthetic court event hearing date
#'
#' @format A date variable.
"hearingdate_synth"

#' Interim Dane data with failure to appear (FTA) as an outcome
#'
#' An interim dataset containing pre-treatment covariates, a binary treatment (Z), an ordinal decision (D),
#' and an outcome variable (Y). The data used for the paper, and made available here, are interim, based on
#' only half of the observations in the study and (for those observations) only half of the study follow-up period.
#' We use them only to illustrate methods, not to draw substantive conclusions.
#'
#' @format A data frame with 1891 rows and 19 variables:
#' \describe{
#'   \item{Z}{binary treatment}
#'   \item{D}{ordinal decision}
#'   \item{Y}{outcome}
#'   \item{Sex}{male or female}
#'   \item{White}{white or non-white}
#'   \item{SexWhite}{the interaction between gender and race}
#'   \item{Age}{age}
#'   \item{PendingChargeAtTimeOfOffense}{binary variable for pending charge (felony, misdemeanor, or both) at the time of offense}
#'   \item{NCorNonViolentMisdemeanorCharge}{binary variable for current non-violent felony charge}
#'   \item{ViolentMisdemeanorCharge}{binary variable for current violent misdemeanor charge}
#'   \item{ViolentFelonyCharge}{binary variable for current violent felony charge}
#'   \item{NonViolentFelonyCharge}{binary variable for current non-violent felony charge}
#'   \item{PriorMisdemeanorConviction}{binary variable for prior conviction of misdemeanor}
#'   \item{PriorFelonyConviction}{binary variable for prior conviction of felony}
#'   \item{PriorViolentConviction}{four-level ordinal variable for prior violent conviction}
#'   \item{PriorSentenceToIncarceration}{binary variable for prior sentence to incarceration}
#'   \item{PriorFTAInPastTwoYears}{three-level ordinal variable for FTAs from past two years}
#'   \item{PriorFTAOlderThanTwoYears}{binary variable for FTAs from over two years ago}
#'   \item{Staff_ReleaseRecommendation}{four-level ordinal variable for the DMF recommendation}
#' }
"FTAdata"

#' Interim Dane data with new criminal activity (NCA) as an outcome
#'
#' An interim dataset containing pre-treatment covariates, a binary treatment (Z), an ordinal decision (D),
#' and an outcome variable (Y). The data used for the paper, and made available here, are interim, based on
#' only half of the observations in the study and (for those observations) only half of the study follow-up period.
#' We use them only to illustrate methods, not to draw substantive conclusions.
#'
#' @format A data frame with 1891 rows and 19 variables:
#' \describe{
#'   \item{Z}{binary treatment}
#'   \item{D}{ordinal decision}
#'   \item{Y}{outcome}
#'   \item{Sex}{male or female}
#'   \item{White}{white or non-white}
#'   \item{SexWhite}{the interaction between gender and race}
#'   \item{Age}{age}
#'   \item{PendingChargeAtTimeOfOffense}{binary variable for pending charge (felony, misdemeanor, or both) at the time of offense}
#'   \item{NCorNonViolentMisdemeanorCharge}{binary variable for current non-violent felony charge}
#'   \item{ViolentMisdemeanorCharge}{binary variable for current violent misdemeanor charge}
#'   \item{ViolentFelonyCharge}{binary variable for current violent felony charge}
#'   \item{NonViolentFelonyCharge}{binary variable for current non-violent felony charge}
#'   \item{PriorMisdemeanorConviction}{binary variable for prior conviction of misdemeanor}
#'   \item{PriorFelonyConviction}{binary variable for prior conviction of felony}
#'   \item{PriorViolentConviction}{four-level ordinal variable for prior violent conviction}
#'   \item{PriorSentenceToIncarceration}{binary variable for prior sentence to incarceration}
#'   \item{PriorFTAInPastTwoYears}{three-level ordinal variable for FTAs from past two years}
#'   \item{PriorFTAOlderThanTwoYears}{binary variable for FTAs from over two years ago}
#'   \item{Staff_ReleaseRecommendation}{four-level ordinal variable for the DMF recommendation}
#' }
"NCAdata"

#' Interim Dane data with new violent criminal activity (NVCA) as an outcome
#'
#' An interim dataset containing pre-treatment covariates, a binary treatment (Z), an ordinal decision (D),
#' and an outcome variable (Y). The data used for the paper, and made available here, are interim, based on
#' only half of the observations in the study and (for those observations) only half of the study follow-up period.
#' We use them only to illustrate methods, not to draw substantive conclusions.
#'
#' @format A data frame with 1891 rows and 19 variables:
#' \describe{
#'   \item{Z}{binary treatment}
#'   \item{D}{ordinal decision}
#'   \item{Y}{outcome}
#'   \item{Sex}{male or female}
#'   \item{White}{white or non-white}
#'   \item{SexWhite}{the interaction between gender and race}
#'   \item{Age}{age}
#'   \item{PendingChargeAtTimeOfOffense}{binary variable for pending charge (felony, misdemeanor, or both) at the time of offense}
#'   \item{NCorNonViolentMisdemeanorCharge}{binary variable for current non-violent felony charge}
#'   \item{ViolentMisdemeanorCharge}{binary variable for current violent misdemeanor charge}
#'   \item{ViolentFelonyCharge}{binary variable for current violent felony charge}
#'   \item{NonViolentFelonyCharge}{binary variable for current non-violent felony charge}
#'   \item{PriorMisdemeanorConviction}{binary variable for prior conviction of misdemeanor}
#'   \item{PriorFelonyConviction}{binary variable for prior conviction of felony}
#'   \item{PriorViolentConviction}{four-level ordinal variable for prior violent conviction}
#'   \item{PriorSentenceToIncarceration}{binary variable for prior sentence to incarceration}
#'   \item{PriorFTAInPastTwoYears}{three-level ordinal variable for FTAs from past two years}
#'   \item{PriorFTAOlderThanTwoYears}{binary variable for FTAs from over two years ago}
#'   \item{Staff_ReleaseRecommendation}{four-level ordinal variable for the DMF recommendation}
#' }
"NVCAdata"

#' Interim Dane PSA data
#'
#' An interim dataset containing a binary treatment (Z), ordinal decision (D), three PSA variables
#' (FTAScore, NCAScore, and NVCAFlag), DMF recommendation, and two pre-treatment covariates (binary
#' indicator for gender; binary indicator for race). The data used for the paper, and made available here, are interim, based on
#' only half of the observations in the study and (for those observations) only half of the study follow-up period.
#' We use them only to illustrate methods, not to draw substantive conclusions.
#'
#' @format A data frame with 1891 rows and 7 variables:
#' \describe{
#'   \item{Z}{binary treatment}
#'   \item{D}{ordinal decision}
#'   \item{FTAScore}{FTA score}
#'   \item{NCAScore}{NCA score}
#'   \item{NVCAFlag}{NVCA flag}
#'   \item{DMF}{DMF recommendation}
#'   \item{Sex}{male or female}
#'   \item{White}{white or non-white}
#' }
"PSAdata"

#' Interim court event hearing date
#'
#' An Interim Dane court event hearing date of Dane data in factor format. The data used for the paper, and made available here, are interim, based on
#' only half of the observations in the study and (for those observations) only half of the study follow-up period.
#' We use them only to illustrate methods, not to draw substantive conclusions.
#'
#' @format A date variable in factor format.
"HearingDate"
