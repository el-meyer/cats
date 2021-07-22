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
#' random <- TRUE
#'
#' rr_comb <- c(0.40, 0.45, 0.50)
#' prob_comb_rr <- c(0.4, 0.4, 0.2)
#' rr_mono <- c(0.20, 0.25, 0.30)
#' prob_mono_rr <- c(0.2, 0.4, 0.4)
#' rr_back <- c(0.20, 0.25, 0.30)
#' prob_back_rr <- c(0.2, 0.4, 0.4)
#' rr_plac <- c(0.10, 0.12, 0.14)
#' prob_plac_rr <- c(0.25, 0.5, 0.25)
#'
#' rr_transform <- list(
#'   function(x) {return(c(0.75*(1 - x), (1-0.75)*(1-x), (1-0.75)*x, 0.75*x))},
#'   function(x) {return(c(0.85*(1 - x), (1-0.85)*(1-x), (1-0.85)*x, 0.85*x))}
#' )
#' prob_rr_transform <- c(0.5, 0.5)
#'
#' cohorts_max <- 4
#' safety_prob <- 0
#' sharing_type <- "all"
#' trial_struc <- "all_plac"
#' sr_drugs_pos <- 4
#' n_int <- 100
#' n_fin <- 200
#' stage_data <- TRUE
#' cohort_random <- 0.05
#' target_rr <- c(0,0,1)
#' cohort_offset <- 0
#' random_type <- "absolute"
#' sr_first_pos <- FALSE
#' missing_prob <- 0.1
#'
#' # Vergleich Combo vs Mono
#' Bayes_Sup1 <- matrix(nrow = 3, ncol = 3)
#' Bayes_Sup1[1,] <- c(0.05, 0.90, 1.00)
#' Bayes_Sup1[2,] <- c(0.05, 0.65, 1.00)
#' Bayes_Sup1[3,] <- c(0.10, 0.50, 1.00)
#' # Vergleich Combo vs Backbone
#' Bayes_Sup2 <- matrix(nrow = 3, ncol = 3)
#' Bayes_Sup2[1,] <- c(0.05, 0.90, 1.00)
#' Bayes_Sup2[2,] <- c(NA, NA, NA)
#' Bayes_Sup2[3,] <- c(NA, NA, NA)
#' # Vergleich Mono vs Placebo
#' Bayes_Sup3 <- matrix(nrow = 3, ncol = 3)
#' Bayes_Sup3[1,] <- c(0.00, 0.90, 1.00)
#' Bayes_Sup3[2,] <- c(NA, NA, NA)
#' Bayes_Sup3[3,] <- c(NA, NA, NA)
#' # Vergleich Back vs Placebo
#' Bayes_Sup4 <- matrix(nrow = 3, ncol = 3)
#' Bayes_Sup4[1,] <- c(0.00, 0.90, 1.00)
#' Bayes_Sup4[2,] <- c(NA, NA, NA)
#' Bayes_Sup4[3,] <- c(NA, NA, NA)
#' Bayes_Sup <- list(list(Bayes_Sup1, Bayes_Sup2, Bayes_Sup3, Bayes_Sup4),
#'              list(Bayes_Sup1, Bayes_Sup2, Bayes_Sup3, Bayes_Sup4))
#'
#' ocs <- trial_ocs(
#' n_int = n_int, n_fin = n_fin, random_type = random_type,
#' rr_comb = rr_comb, rr_mono = rr_mono, rr_back = rr_back, rr_plac = rr_plac,
#' rr_transform = rr_transform, random = random, prob_comb_rr = prob_comb_rr,
#' prob_mono_rr = prob_mono_rr, prob_back_rr = prob_back_rr, prob_plac_rr = prob_plac_rr,
#' stage_data = stage_data, cohort_random = cohort_random, cohorts_max = cohorts_max,
#' sr_drugs_pos = sr_drugs_pos, target_rr = target_rr, sharing_type = sharing_type,
#' sr_first_pos = sr_first_pos, safety_prob = safety_prob, Bayes_Sup = Bayes_Sup,
#' prob_rr_transform = prob_rr_transform, cohort_offset = cohort_offset,
#' trial_struc = trial_struc, missing_prob = missing_prob,
#' iter = 150, coresnum = 1, save = FALSE, ret_list = TRUE, plot_ocs = TRUE
#' )
#'
#' ocs[[3]]
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
    trial_results <- foreach::foreach(i = 1:iter, .packages = "CohortPlat", .export = export) %dopar% {
      # first call program function
      trial_res <- do.call(simulate_trial, arguments)
      # Now save individual trial results
      trial_res
    }
    # end parallel
    doParallel::stopImplicitCluster()
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
    Avg_Pat_First_Suc          = mean(sapply(trial_results, function(x) x$Trial_Overview$Total_N_First_Suc), na.rm = TRUE),
    Avg_Perc_Pat_Sup_Plac_Th   = mean(sapply(trial_results, function(x) x$Trial_Overview$Perc_N_Sup_Plac_Th)),
    Avg_Perc_Pat_Sup_Plac_Real = mean(sapply(trial_results, function(x) x$Trial_Overview$Perc_N_Sup_Plac_Real)),
    Avg_Pat_Comb               = mean(sapply(trial_results, function(x) x$Trial_Overview$Total_N_Comb)),
    Avg_Pat_Mono               = mean(sapply(trial_results, function(x) x$Trial_Overview$Total_N_Mono)),
    Avg_Pat_Back               = mean(sapply(trial_results, function(x) x$Trial_Overview$Total_N_Back)),
    Avg_Pat_Plac               = mean(sapply(trial_results, function(x) x$Trial_Overview$Total_N_Plac)),
    Avg_Pat_Plac_First_Suc     = mean(sapply(trial_results, function(x) x$Trial_Overview$Total_N_Plac_First_Suc), na.rm = TRUE),
    Avg_Pat_Plac_Pool          = mean(sapply(trial_results, function(x) x$Trial_Overview$Total_N_Plac_Pool), na.rm = TRUE),
    Avg_RR_Comb                = mean(unlist(sapply(trial_results, function(x) x$Trial_Overview$RR_Comb)), na.rm = TRUE),
    Avg_RR_Mono                = mean(unlist(sapply(trial_results, function(x) x$Trial_Overview$RR_Mono)), na.rm = TRUE),
    Avg_RR_Back                = mean(unlist(sapply(trial_results, function(x) x$Trial_Overview$RR_Back)), na.rm = TRUE),
    Avg_RR_Plac                = mean(unlist(sapply(trial_results, function(x) x$Trial_Overview$RR_Plac)), na.rm = TRUE),
    SD_RR_Comb                 = stats::sd(unlist(sapply(trial_results, function(x) x$Trial_Overview$RR_Comb)), na.rm = TRUE),
    SD_RR_Mono                 = stats::sd(unlist(sapply(trial_results, function(x) x$Trial_Overview$RR_Mono)), na.rm = TRUE),
    SD_RR_Back                 = stats::sd(unlist(sapply(trial_results, function(x) x$Trial_Overview$RR_Back)), na.rm = TRUE),
    SD_RR_Plac                 = stats::sd(unlist(sapply(trial_results, function(x) x$Trial_Overview$RR_Plac)), na.rm = TRUE),
    Avg_Suc_Hist               = mean(sapply(trial_results, function(x) x$Trial_Overview$Successes_Hist)),
    Avg_Suc_Hist_Comb          = mean(sapply(trial_results, function(x) x$Trial_Overview$Successes_Hist_Comb)),
    Avg_Suc_Hist_Mono          = mean(sapply(trial_results, function(x) x$Trial_Overview$Successes_Hist_Mono)),
    Avg_Suc_Hist_Back          = mean(sapply(trial_results, function(x) x$Trial_Overview$Successes_Hist_Back)),
    Avg_Suc_Hist_Plac          = mean(sapply(trial_results, function(x) x$Trial_Overview$Successes_Hist_Plac)),
    Avg_Suc_Bio                = mean(sapply(trial_results, function(x) x$Trial_Overview$Successes_Bio)),
    Avg_Suc_Bio_Comb           = mean(sapply(trial_results, function(x) x$Trial_Overview$Successes_Bio_Comb)),
    Avg_Suc_Bio_Mono           = mean(sapply(trial_results, function(x) x$Trial_Overview$Successes_Bio_Mono)),
    Avg_Suc_Bio_Back           = mean(sapply(trial_results, function(x) x$Trial_Overview$Successes_Bio_Back)),
    Avg_Suc_Bio_Plac           = mean(sapply(trial_results, function(x) x$Trial_Overview$Successes_Bio_Plac)),
    Avg_Cohorts                = mean(sapply(trial_results, function(x) x$Trial_Overview$N_Cohorts)),
    Avg_Cohorts_First_Suc      = mean(sapply(trial_results, function(x) x$Trial_Overview$N_Cohorts_First_Suc), na.rm = TRUE),
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
