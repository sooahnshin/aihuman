## code to prepare `nuis_func` and `nuis_func_ai`

library(tidyr)
library(dplyr)
library(aihuman)

# randomized PSA provision (0: none, 1: provided)
Z <- aihuman::NCAdata$Z
# judge's release decision (0: signature, 1: cash)
D <- if_else(aihuman::NCAdata$D == 0, 0, 1)
# dichotomized pretrial public safety assessment scores (0: signature, 1: cash)
A <- aihuman::PSAdata$DMF
# new criminal activity (0: no, 1: yes)
Y <- aihuman::NCAdata$Y
# covariates
cov_mat <- aihuman::NCAdata |> 
  select(-c(Y, D, Z)) |> 
  bind_cols(FTAScore = aihuman::PSAdata$FTAScore, 
            NCAScore = aihuman::PSAdata$NCAScore, 
            NVCAFlag = aihuman::PSAdata$NVCAFlag) |>
  as.matrix()

set.seed(12345)
nuis_func <- compute_nuisance_functions(Y = Y, D = D, Z = Z, V = cov_mat, shrinkage = 0.01, n.trees = 1000)
nuis_func_ai <- compute_nuisance_functions_ai(Y = Y, D = D, Z = Z, A = A, V = cov_mat, shrinkage = 0.01, n.trees = 1000)

usethis::use_data(nuis_func, nuis_func_ai, internal = TRUE)
