#' Calculates the operating characteristics of the cohort trial
#'
#' Given the trial specific design parameters, performs a number
#' of simulations of the trial and saves the result in an Excel file
#'
#' @param iter         Number of program simulations that should be performed
#'
#' @param coresnum     How many cores should be used for parallel computing
#'
#' @param save         Indicator whether simulation results should be saved in an Excel file
#'
#' @param path         Path to which simulation results will be saved; if NULL, then save to current path
#'
#' @param ret_list     Indicator whether function should return list of results
#'
#' @param ret_trials   Indicator whether individual trial results should be saved as well
#'
#' @param filename     Filename of saved Excel file with results; if NULL, then name will contain design parameters
#'
#' @param plot_ocs     Should OCs stability plots be drawn?
#'
#' @param export       Should any other variables be exported to the parallel tasks?
#'
#' @param ...          All other design parameters for chosen program
#'
#' @return List containing: Responses and patients on experimental and control arm, total treatment successes and failures and final p-value
#'
#' @examples
#'
# random <- TRUE
#
# rr_comb1 <- 0.25
# prob_comb1_rr <- 1
# rr_comb2 <- 0.20
# prob_comb2_rr <- 1
# rr_plac1 <- 0.10
# prob_plac1_rr <- 1
# rr_plac2 <- 0.10
# prob_plac2_rr <- 1
#
# correlation <- 0.8
#
# cohorts_start <- 2
# cohorts_max <- 5
# safety_prob <- 0
# sharing_type <- "concurrent"
# sr_drugs_pos <- 5
# sr_first_pos <- FALSE
# n_fin <- 100
# stage_data <- TRUE
# cohort_random <- 0.01
# cohort_offset <- 0
# cohorts_sim <- Inf
# random_type <- "absolute"
# missing_prob <- 0.2
# cohort_fixed <- 5
# hist_lag <- 48
# analysis_times <- c(0.5, 0.75, 1)
# accrual_type <- "fixed"
# accrual_param <- 9
# time_trend <- 0.001
# composite <- "or"
#'
#' # Comparison IA1
#' Bayes_Sup11 <- matrix(nrow = 2, ncol = 2)
#' Bayes_Sup11[1,] <- c(0.00, 0.95)
#' Bayes_Sup11[2,] <- c(0.10, 0.80)
#' # Comparison IA2
#' Bayes_Sup12 <- matrix(nrow = 2, ncol = 2)
#' Bayes_Sup12[1,] <- c(0.00, 0.95)
#' Bayes_Sup12[2,] <- c(0.10, 0.80)
#' # Comparison IA3
#' Bayes_Sup13 <- matrix(nrow = 2, ncol = 2)
#' Bayes_Sup13[1,] <- c(0.00, 0.95)
#' Bayes_Sup13[2,] <- c(0.10, 0.80)
#'
#' Bayes_Sup1 <- Bayes_Sup2 <- list(list(Bayes_Sup11), list(Bayes_Sup12), list(Bayes_Sup13))
#'
# ocs <- trial_ocs(
# n_fin = n_fin, random_type = random_type, composite = composite,
# rr_comb1 = rr_comb1, rr_comb2 = rr_comb2, rr_plac1 = rr_plac1, rr_plac2 = rr_plac2,
# random = random, prob_comb1_rr = prob_comb1_rr, prob_comb2_rr = prob_comb2_rr,
# prob_plac1_rr = prob_plac1_rr, prob_plac2_rr = prob_plac2_rr,
# stage_data = stage_data, cohort_random = cohort_random, cohorts_max = cohorts_max,
# sr_drugs_pos = sr_drugs_pos, sharing_type = sharing_type, correlation = correlation,
# safety_prob = safety_prob, Bayes_Sup1 = Bayes_Sup1, Bayes_Sup2 = Bayes_Sup2,
# cohort_offset = cohort_offset, sr_first_pos = sr_first_pos,
# missing_prob = missing_prob, cohort_fixed = cohort_fixed, accrual_type = accrual_type,
# accrual_param = accrual_param, hist_lag = hist_lag, analysis_times = analysis_times,
# time_trend = time_trend, cohorts_start = cohorts_start, cohorts_sim = cohorts_sim,
# iter = 2, coresnum = 1, save = FALSE, ret_list = TRUE, plot_ocs = TRUE
# )
#'
#' @export
trial_ocs <- function(iter, coresnum = 1, save = FALSE, path = NULL, filename = NULL, ret_list = FALSE,
                         ret_trials = FALSE, plot_ocs = FALSE, export = NULL, ...) {

  ##### Initialize variables #####

  # Since R CMD check allows only for 2 cores, set this
  if (coresnum > 1) {
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(chk) && chk == "TRUE") {
      coresnum <- 2
    }

    # Prepare for parallel computing
    cl <- parallel::makePSOCKcluster(coresnum)
    doParallel::registerDoParallel(cl)

    "%dopar%" <- foreach::"%dopar%"

    arguments <- list(...) # gather additional program arguments

    ##### Run parallel simulations #####

    # run in parallel
    trial_results <- foreach::foreach(i = 1:iter, .packages = "cats", .export = export) %dopar% {
      # first call program function
      trial_res <- do.call(cats::simulate_trial, arguments)
      # Now save individual trial results
      trial_res
    }
    # end parallel
    # closeAllConnections()
    parallel::stopCluster(cl)

  } else {

    arguments <- list(...) # gather additional program arguments

    ##### Run parallel simulations #####

    # run without parallel
    trial_results <- list()

    for (i in 1:iter) {
      # first call program function
      trial_res <- do.call(simulate_trial, arguments)
      # Now save individual trial results
      trial_results <- c(trial_results, list(trial_res))
    }

  }



  ##### Compute OCs #####

  # Return OCs
  ret1 <- list(
    Avg_Pat                    = mean(sapply(trial_results, function(x) x$Trial_Overview$Total_N)),
    Avg_Pat_First_Suc          = mean(sapply(trial_results, function(x) x$Trial_Overview$Pat_First_Suc), na.rm = TRUE),
    Avg_Cohorts                = mean(sapply(trial_results, function(x) x$Trial_Overview$N_Cohorts)),
    Avg_Cohort_N               = mean(sapply(trial_results, function(x) x$Trial_Overview$Final_N_Cohort_Trial)),
    Avg_Time                   = mean(sapply(trial_results, function(x) x$Trial_Overview$Total_Time)),
    Avg_Time_First_Suc         = mean(sapply(trial_results, function(x) x$Trial_Overview$Time_First_Suc), na.rm = TRUE),
    Avg_Unenrolled_Pat         = mean(sapply(trial_results, function(x) x$Trial_Overview$Unenrolled_Pats), na.rm = TRUE),

    Avg_Intx1_Go                = sum(sapply(trial_results, function(x) x$Trial_Overview$Intx1_GO)) /
                                  sum(sapply(trial_results, function(x) x$Trial_Overview$N_Cohorts)),
    Avg_Intx2_Go                = sum(sapply(trial_results, function(x) x$Trial_Overview$Intx2_GO)) /
      sum(sapply(trial_results, function(x) x$Trial_Overview$N_Cohorts)),
    Avg_Intx3_Go                = sum(sapply(trial_results, function(x) x$Trial_Overview$Intx3_GO)) /
      sum(sapply(trial_results, function(x) x$Trial_Overview$N_Cohorts)),

    Avg_Intx1_Stop              = sum(sapply(trial_results, function(x) x$Trial_Overview$Intx1_STOP)) /
      sum(sapply(trial_results, function(x) x$Trial_Overview$N_Cohorts)),
    Avg_Intx2_Stop              = sum(sapply(trial_results, function(x) x$Trial_Overview$Intx2_STOP)) /
      sum(sapply(trial_results, function(x) x$Trial_Overview$N_Cohorts)),
    Avg_Intx3_Stop              = sum(sapply(trial_results, function(x) x$Trial_Overview$Intx3_STOP)) /
      sum(sapply(trial_results, function(x) x$Trial_Overview$N_Cohorts)),

    Avg_Rand_Stop              = sum(sapply(trial_results, function(x) x$Trial_Overview$Safety_STOP)) /
      sum(sapply(trial_results, function(x) x$Trial_Overview$N_Cohorts)),

    Avg_TP                     = mean(sapply(trial_results, function(x) x$Trial_Overview$TP)),
    Avg_FP                     = mean(sapply(trial_results, function(x) x$Trial_Overview$FP)),
    Avg_TN                     = mean(sapply(trial_results, function(x) x$Trial_Overview$TN)),
    Avg_FN                     = mean(sapply(trial_results, function(x) x$Trial_Overview$FN)),
    Avg_any_P                  = mean(sapply(trial_results, function(x) x$Trial_Overview$any_P)),
    Dist_FWER                  = sapply(trial_results, function(x) x$Trial_Overview$FP > 0),
    Dist_FDR                   = cumsum(sapply(trial_results, function(x) x$Trial_Overview$FP)) /
                                 cumsum(sapply(trial_results, function(x) x$Trial_Overview$TP) + sapply(trial_results, function(x) x$Trial_Overview$FP)),
    Dist_Disj_Power            = sapply(trial_results, function(x) x$Trial_Overview$TP > 0),
    Dist_PTP                   = cumsum(sapply(trial_results, function(x) x$Trial_Overview$TP)) /
                                 cumsum((sapply(trial_results, function(x) x$Trial_Overview$TP) + sapply(trial_results, function(x) x$Trial_Overview$FN))),
    Dist_PTT1ER                = cumsum(sapply(trial_results, function(x) x$Trial_Overview$FP)) /
                                 cumsum((sapply(trial_results, function(x) x$Trial_Overview$FP) + sapply(trial_results, function(x) x$Trial_Overview$TN))),
    FDR                        = sum(sapply(trial_results, function(x) x$Trial_Overview$FP)) /
                                 sum(sapply(trial_results, function(x) x$Trial_Overview$TP) + sapply(trial_results, function(x) x$Trial_Overview$FP)),
    PTP                        = sum(sapply(trial_results, function(x) x$Trial_Overview$TP)) /
                                 sum((sapply(trial_results, function(x) x$Trial_Overview$TP) + sapply(trial_results, function(x) x$Trial_Overview$FN))),
    PTT1ER                     = sum(sapply(trial_results, function(x) x$Trial_Overview$FP)) /
                                 sum((sapply(trial_results, function(x) x$Trial_Overview$FP) + sapply(trial_results, function(x) x$Trial_Overview$TN)))
  )

  any_H0 <- sapply(trial_results, function(x) x$Trial_Overview$FP > 0) | sapply(trial_results, function(x) x$Trial_Overview$TN > 0)
  any_H1 <- sapply(trial_results, function(x) x$Trial_Overview$TP > 0) | sapply(trial_results, function(x) x$Trial_Overview$FN > 0)

  # Get "classical" operating characteristics
  ret1$FWER <- mean(as.numeric(ret1$Dist_FWER)[any_H0])
  ret1$Disj_Power <- mean(as.numeric(ret1$Dist_Disj_Power)[any_H1])

  # Get "Bayesian Average" operating characteristics
  ret1$FWER_BA <- mean(ret1$Dist_FWER)
  ret1$Disj_Power_BA <- mean(ret1$Dist_Disj_Power)

  # Get all function arguments to display later
  arguments_full <- list(iter = iter, coresnum = coresnum, ...)

  ##### Should Simulation be plotted? ######

  if (plot_ocs) {
    "%>%" <- dplyr::"%>%"

    # Prepare plot

    fwer_classical <- ret1$Dist_FWER
    fwer_classical[!any_H0] <- NA
    fwer_classical_dist <- fwer_classical
    fwer_classical_dist[!is.na(fwer_classical)] <- dplyr::cummean(fwer_classical_dist[!is.na(fwer_classical)])
    fwer_classical_dist <- zoo::na.locf(fwer_classical_dist)
    if (length(fwer_classical_dist) != length(ret1$Dist_FWER)) {
      fwer_classical_dist <- c(rep(0, length(ret1$Dist_FWER) - length(fwer_classical_dist)), fwer_classical_dist)
    }

    disj_power_classical <- ret1$Dist_Disj_Power
    disj_power_classical[!any_H1] <- NA
    disj_power_classical_dist <- disj_power_classical
    disj_power_classical_dist[!is.na(disj_power_classical)] <- dplyr::cummean(disj_power_classical_dist[!is.na(disj_power_classical)])
    disj_power_classical_dist <- zoo::na.locf(disj_power_classical_dist)
    if (length(disj_power_classical_dist) != length(ret1$Dist_FWER)) {
      disj_power_classical_dist <- c(rep(0, length(ret1$Dist_FWER) - length(disj_power_classical_dist)), disj_power_classical_dist)
    }

    d1 <- dplyr::tibble(
      FWER_BA = dplyr::cummean(ret1$Dist_FWER),
      Disj_Power_BA = dplyr::cummean(ret1$Dist_Disj_Power),
      FWER_CD = fwer_classical_dist,
      Disj_Power_CD = disj_power_classical_dist,
      FDR = ret1$Dist_FDR,
      PTP = ret1$Dist_PTP,
      PTT1ER = ret1$Dist_PTT1ER
      ) %>%
      tidyr::gather(
        key = "Error_Rate", value = "Prob",
        FWER_CD, FWER_BA, FDR, Disj_Power_CD, Disj_Power_BA, PTP, PTT1ER,
        factor_key = TRUE
      ) %>%
      dplyr::mutate(
        Simulation = rep(1:iter, 7)
      )

    sim_plot <-
      plotly::ggplotly(
        ggplot2::ggplot(d1, ggplot2::aes(x = Simulation, y = Prob, color = Error_Rate)) +
          ggplot2::geom_line() +
          ggplot2::theme_minimal()
      )
    sim_plot$x$layout$annotations[[1]]$text <- ""
  }

  ##### Save as Excel and RData #####

  if (save) {
    # If results should be saved, save to Excel and slightly recode ret1 to be a matrix
    # If path is supplied, results will be saved in folder with program name at this path; if folder with program names does not exist, create one
    # If path is not supplied, go to current WD, create folder "temp" and proceed analogously

    # Get the above results (which are a list) and convert to 1xk vector for better display in excel file
    ret2 <- t(as.matrix(ret1))
    # Create return object which includes a sheet (==list element) with all the design parameters,
    # the program OCs and all the program simulations results

    arguments_full2 <- unlist(arguments_full)
    arguments_full2[which(sapply(arguments_full2, function(x) is.function(x)))] <-
      as.character(arguments_full2[which(sapply(arguments_full2, function(x) is.function(x)))])
    ret <- list(t(arguments_full2), ret2)

    # Additionall check whether Unix (Mac, Linux) or not to account for minor differences
    if (.Platform$OS.type == "unix") {
      if (is.null(path)) {
        path0 <- getwd()
        ifelse(!dir.exists(file.path(path0, "tempsim/")), dir.create(file.path(path0, "temp/")), FALSE)
        path <- file.path(path0, "tempsim/")
      }

      file.savepath <- paste0(path, filename, ".xlsx")

      # Write to xlsx
      openxlsx::write.xlsx(ret, file = file.savepath)
      # Save as RData
      file.savepath.rdata <- paste0(path, filename, ".RData")
      results <- c(list(arguments_full), list(ret1), trial_results)
      save(results, file = file.savepath.rdata)

    } else {
      if (is.null(path)) {
        path0 <- getwd()
        ifelse(!dir.exists(file.path(path0, "tempsim")), dir.create(file.path(path0, "tempsim")), FALSE)
        path <- file.path(path0, "tempsim")
      }

      file.savepath <- paste0(path, "/", filename, ".xlsx")

      # Write to xlsx
      openxlsx::write.xlsx(ret, file = file.savepath)
      # Save as RData
      file.savepath.rdata <- paste0(path, "/", filename, ".RData")
      results <- c(list(arguments_full), list(ret1), trial_results)
      save(results, file = file.savepath.rdata)
    }
  }

  ##### Return Values #####

  # If result should be also returned as a list, do so
  if (ret_list) {
    if (plot_ocs) {
      ret <- c(list(arguments_full), list(ret1), list(sim_plot))
    } else {
      ret <- c(list(arguments_full), list(ret1))
    }
    # if also individual trial data should be returned, add that to list
    if (ret_trials) {
      if (plot_ocs) {
        ret <- c(list(arguments_full), list(ret1), trial_results, list(sim_plot))
      } else {
        ret <- c(list(arguments_full), list(ret1), trial_results)
      }
    }
    return(ret)
  }
}
