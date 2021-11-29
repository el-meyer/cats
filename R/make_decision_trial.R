#' Checks whether decision criteria are met and updates trial results accordingly.
#'
#' Given a res_list object, checks the supplied decision criteria and saves the results in the res_list file.
#'
#' @param res_list         List item containing individual cohort trial results so far in a format used by the
#'                         other functions in this package
#'
#' @param which_cohort     Current cohort that should be evaluated
#'
#' @param w                If dynamic borrowing, what is the prior choice for w. Default is 0.5.
#'
#' @param beta_prior       Prior parameter for all Beta Distributions. Default is 0.5.
#'
#' @param Bayes_Sup1       List of matrices with rows corresponding to number of multiple Bayesian posterior two-arm combination criteria for superiority of histology endpoint 1
#'
#' @param Bayes_Sup2       List of matrices with rows corresponding to number of multiple Bayesian posterior two-arm combination criteria for superiority of histology endpoint 2
#'
#' @param Bayes_Fut1       List of matrices with rows corresponding to number of multiple Bayesian posterior two-arm combination criteria for futility of histology endpoint 1
#'
#' @param Bayes_Fut2       List of matrices with rows corresponding to number of multiple Bayesian posterior two-arm combination criteria for futility of histology endpoint 2
#'
#' @param analysis_number  1st, second or third analysis?
#'
#' @param analysis_time    Platform Time of Analysis
#'
#' @param dataset          Dataset to be used for analysis
#'
#' @param hist_miss        Whether or not to exclude missing histology data
#'
#' @param hist_lag         Histology Lag
#'
#' @param sharing_type     Type of Data Sharing to perform
#'
#' @param endpoint_number  Should histology endpoint 1 or 2 be evaluated?
#'
#' @param ...              Further arguments inherited from simulate_trial
#'
#' @return List containing original res_list and results of decision rules
#'
#' @examples
#'
#' # Example 1
#'
#' # Initialize empty data frame
#' cols <- c("PatID", "ArrivalTime", "Cohort", "Arm", "RespHist1", "RespHist2", "HistMissing")
#' df <- matrix(nrow = 100, ncol = length(cols))
#' colnames(df) <- cols
#' df <- as.data.frame(df)
#' df$PatID <- 1:100
#' df$ArrivalTime <- sort(runif(100, min = 0, max = 5))
#' df$Cohort <- sample(1:2, 100, replace = TRUE)
#' df$Arm <- sample(c("Combo", "Plac"), 100, replace = TRUE)
#' df$RespHist1 <- sample(0:1, 100, replace = TRUE)
#' df$RespHist2 <- sample(0:1, 100, replace = TRUE)
#' df$HistMissing <- sample(0:1, 100, replace = TRUE, prob = c(0.95, 0.05))
#'
#' # Initialize res_list Object
#'
#' res_list <-
#' rep(
#'   list(
#'     list(
#'       Meta = list(
#'         decision = rep("none", 3),
#'         decision_hist1 = rep("none", 3),
#'         decision_hist2 = rep("none", 3),
#'         start_n = 0,
#'         start_time = 0,
#'         pat_enrolled = 0
#'       ),
#'       Arms = rep(
#'         list(
#'           list(
#'             rr = NULL,
#'             hist_observed = 0
#'           )
#'         ),
#'         2
#'       )
#'     )
#'   ),
#'   2
#' )
#'
#' arm_names <- c("Comb", "Plac")
#'
#' for (i in 1:2) {
#'   names(res_list)[i] <- paste0("Cohort", i)
#'   names(res_list[[i]]$Arms) <- arm_names
#'
#'   res_list[[i]]$Arms$Comb$rr <- matrix(c(0.2, 0.2), ncol = 2)
#'   res_list[[i]]$Arms$Plac$rr <- matrix(c(0.1, 0.1), ncol = 2)
#' }
#'
#' sharing_type <- "all"
#' analysis_number <- 3
#' which_cohort <- 1
#' endpoint_number <- 2
#' hist_lag <- 1
#' analysis_time <- 6
#'
#' # Comparison IA1
#' Bayes_Sup11 <- matrix(nrow = 2, ncol = 2)
#' Bayes_Sup11[1,] <- c(0.00, 0.95)
#' Bayes_Sup11[2,] <- c(0.10, 0.80)
#' # Comparison IA2
#' Bayes_Sup12 <- matrix(nrow = 2, ncol = 2)
#' Bayes_Sup12[1,] <- c(0.00, 0.95)
#' Bayes_Sup12[2,] <- c(NA, NA)
#' # Comparison IA3
#' Bayes_Sup13 <- matrix(nrow = 2, ncol = 2)
#' Bayes_Sup13[1,] <- c(0.00, 0.95)
#' Bayes_Sup13[2,] <- c(0.10, 0.80)
#'
#' Bayes_Sup1 <- Bayes_Sup2 <- list(list(Bayes_Sup11), list(Bayes_Sup12), list(Bayes_Sup13))
#'
#'
#' # DO NOT RUN
#' res_list2 <-
#' make_decision_trial(
#' res_list = res_list, which_cohort = which_cohort,
#' analysis_number = analysis_number, endpoint_number = endpoint_number,
#' Bayes_Sup1 = Bayes_Sup1, Bayes_Sup2 = Bayes_Sup2,
#' dataset = df, analysis_time = analysis_time, hist_lag = hist_lag,
#' sharing_type = sharing_type
#' )
#'
#' @export
make_decision_trial <- function(res_list, which_cohort, Bayes_Sup1 = NULL, Bayes_Fut1 = NULL,
                                Bayes_Sup2 = NULL, Bayes_Fut2 = NULL,
                                w = 0.5, analysis_number, beta_prior = 0.5, hist_lag,
                                endpoint_number, analysis_time, dataset,
                                hist_miss = TRUE, sharing_type, ...) {

  ########## Parameter Explanation and Variable and Function Initiation ##########

  # Bayes_Sup = Matrix with rows corresponding to number of multiple Bayesian two-arm posterior combination criteria for superiority
  # First element is always the delta, second element is the GO-threshold
  # The rule is: If P(P_e > P_c + delta) > GO-threshold, then GO
  # The second rule is: If promising_threshold < P(P_e > P_c + delta) < GO-threshold, then declare promising

  # Bayes_Fut = Matrix with rows corresponding to number of multiple Bayesian two-arm posterior combination criteria for futility
  # First element is always the delta and the second element is the STOP-threshold
  # The rule is: If P(P_e > P_c + delta) < STOP-threshold, then STOP


  # Helper function posterior probability of one Beta being by margin delta greater than other Beta
  post_prob_bin <- function(n_exp, n_contr, resp_exp, resp_contr, delta,
                            a0_exp, b0_exp, a0_contr, b0_contr) {

    # in notation of diploma thesis, this calculates the probability P(P_e >= P_c + delta_sup)
    prob_sup <- stats::integrate(function(y) {
      stats::dbeta(y, a0_exp + resp_exp, b0_exp - resp_exp + n_exp) *
        stats::pbeta(y - delta, a0_contr + resp_contr, b0_contr - resp_contr + n_contr)
    }, delta, 1)$value

    # return posterior probability
    return(prob_sup)
  }

  # Create list for whether futility/superiority decisions were made for every single rule for all three comparisons
  # All list elements will be matrices, with rows corresponding to multiple decision rules and columns corresponding to the
  # three comparisons (j=1: Combo vs Mono, j=2: Combo vs Back, j=3: Mono vs Placebo, j=4: Back vs Placebo)
  superiority_reached <- list()
  futility_reached <- list()

  sup_reached <- fut_reached <- NULL

  # Get responders and sample sizes
  resp_hist_tot <- n_tot <- NULL

  # Extract data column

  if (endpoint_number == 1) {
    dat_vec <- dataset$RespHist1
  } else if (endpoint_number == 2) {
    dat_vec <- dataset$RespHist2
  }

  # Combo

    if (hist_miss) {
      ind <-
        which(
          dataset$Cohort == which_cohort &
            dataset$Arm == "Comb" &
            dataset$ArrivalTime < analysis_time - hist_lag &
            !dataset$HistMissing
        )
    } else {
      ind <-
        which(
          dataset$Cohort == which_cohort &
            dataset$Arm == "Comb" &
            dataset$ArrivalTime < analysis_time - hist_lag
        )
    }

    resp_hist <- dat_vec[ind]

    resp_hist_tot[1] <- sum(resp_hist)
    n_tot[1] <- length(resp_hist)


  # Plac

    # if sharing all data, time frame to take under consideration is simply all data from all cohorts
    if (sharing_type == "all") {

      if (hist_miss) {
        ind <-
          which(
              dataset$Arm == "Plac" &
              dataset$ArrivalTime < analysis_time - hist_lag &
              !dataset$HistMissing
          )
      } else {
        ind <-
          which(
              dataset$Arm == "Plac" &
              dataset$ArrivalTime < analysis_time - hist_lag
          )
      }

      resp_hist <- dat_vec[ind]

      resp_hist_tot[2] <- sum(resp_hist)
      n_tot[2] <- length(resp_hist)
    }

    # if sharing no data, take only data from current Arm
    if (sharing_type == "cohort") {

      if (hist_miss) {
        ind <-
          which(
            dataset$Cohort == which_cohort &
              dataset$Arm == "Plac" &
              dataset$ArrivalTime < analysis_time - hist_lag &
              !dataset$HistMissing
          )
      } else {
        ind <-
          which(
            dataset$Cohort == which_cohort &
              dataset$Arm == "Plac" &
              dataset$ArrivalTime < analysis_time - hist_lag
          )
      }

      resp_hist <- dat_vec[ind]

      resp_hist_tot[2] <- sum(resp_hist)
      n_tot[2] <- length(resp_hist)
    }

    # if sharing concurrent data, take only data from certain time period
    if (sharing_type == "concurrent") {

      if (hist_miss) {
        ind <-
          which(
              dataset$Arm == "Plac" &
              dataset$ArrivalTime < analysis_time - hist_lag &
              dataset$ArrivalTime > res_list[[which_cohort]]$Meta$start_time &
              !dataset$HistMissing
          )
      } else {
        ind <-
          which(
              dataset$Arm == "Plac" &
              dataset$ArrivalTime < analysis_time - hist_lag &
              dataset$ArrivalTime > res_list[[which_cohort]]$Meta$start_time
          )
      }

      resp_hist <- dat_vec[ind]

      resp_hist_tot[2] <- sum(resp_hist)
      n_tot[2] <- length(resp_hist)
    }

    # if dynamic borrowing, create two datasets and combine
    if (sharing_type == "dynamic") {

      if (hist_miss) {
        ind_c <-
          which(
            dataset$Cohort == which_cohort &
              dataset$Arm == "Plac" &
              dataset$ArrivalTime < analysis_time - hist_lag &
              !dataset$HistMissing
          )
      } else {
        ind_c <-
          which(
            dataset$Cohort == which_cohort &
              dataset$Arm == "Plac" &
              dataset$ArrivalTime < analysis_time - hist_lag
          )
      }

      if (hist_miss) {
        ind_e <-
          which(
            dataset$Cohort != which_cohort &
              dataset$Arm == "Plac" &
              dataset$ArrivalTime < analysis_time - hist_lag &
              !dataset$HistMissing
          )
      } else {
        ind_e <-
          which(
            dataset$Cohort != which_cohort &
              dataset$Arm == "Plac" &
              dataset$ArrivalTime < analysis_time - hist_lag
          )
      }

      resp_hist_cohort <- dat_vec[ind_c]
      resp_hist_external <- dat_vec[ind_e]

      suc_hist_c <- sum(resp_hist_cohort)
      N_c <- length(resp_hist_cohort)

      # compute historical responders
      suc_hist_h <- sum(resp_hist_external)
      N_h <- length(resp_hist_external)

      w1_hist <- w * beta(suc_hist_c + suc_hist_h + beta_prior, N_h + N_c - suc_hist_c - suc_hist_h + beta_prior) /
        beta(suc_hist_h + beta_prior, N_h - suc_hist_h + beta_prior)
      w2_hist <- (1 - w) * beta(suc_hist_c + beta_prior, N_c - suc_hist_c + beta_prior) / beta(beta_prior, beta_prior)
      w1n_hist <- w1_hist / (w1_hist + w2_hist)
      w2n_hist <- w2_hist / (w1_hist + w2_hist)

      m_a_hist <- w1n_hist * (suc_hist_c + suc_hist_h + beta_prior) + w2n_hist * (suc_hist_c + beta_prior)
      m_b_hist <- w1n_hist * (N_h + N_c - suc_hist_c - suc_hist_h + beta_prior) + w2n_hist * (N_c - suc_hist_c + beta_prior)

      n_tot[2] <- round(m_a_hist + m_b_hist)
      resp_hist_tot[2] <- round(m_a_hist)
    }


  resp_tot <- resp_hist_tot

  ########## Bayesian Two-Arm Superiority Criteria ###############

  if (!is.null(Bayes_Sup1)) {

    # Get the decision rule for interim/final
    if (analysis_number == 1) {
      Bayes_Sup1 <- Bayes_Sup1[[1]]
      Bayes_Sup2 <- Bayes_Sup2[[1]]
    } else if (analysis_number == 2) {
      Bayes_Sup1 <- Bayes_Sup1[[2]]
      Bayes_Sup2 <- Bayes_Sup2[[2]]
    } else if (analysis_number == 3) {
      Bayes_Sup1 <- Bayes_Sup1[[3]]
      Bayes_Sup2 <- Bayes_Sup2[[3]]
    }

    if (endpoint_number == 1) {
      Bayes_Sup <- Bayes_Sup1
    } else if (endpoint_number == 2) {
      Bayes_Sup <- Bayes_Sup2
    }


    # Firstly initiate as many vectors as will maximum be needed (which is determined by the most complex decision rule)
    for (k in 1:max(sapply(Bayes_Sup, function(x) nrow(x)))) {
      assign(paste0("prob_sup_vec", k), rep(NA, 1))
    }

    assign(paste0("sup_reached"), matrix(NA, nrow = max(sapply(Bayes_Sup, function(x) nrow(x))), ncol = 1))

    # For each of the comparisons separately, evaluate the applying decision rules
    # j=1: Combo vs Mono, j=2: Combo vs Back, j=3: Mono vs Placebo

    for (i in 1:max(sapply(Bayes_Sup, function(x) nrow(x)))) {
      for (j in 1:1) {
        # Get which type of comparison should be conducted
        c1 <- 1
        c2 <- 2

        # Evaluate decision rule, but if non-existant, give NA
        if (is.na(Bayes_Sup[[j]][i,1])) {
          eval(parse(text = paste(paste0("prob_sup_vec", i, "[", j, "]"), "<- NA")))
        } else {
          eval(parse(text = paste(paste0("prob_sup_vec", i, "[", j, "]"), "<-",
                                  post_prob_bin(
                                    n_tot[c1], n_tot[c2],
                                    resp_tot[c1], resp_tot[c2],
                                    Bayes_Sup[[j]][i,1], beta_prior, beta_prior, beta_prior, beta_prior))))
        }

        # Check whether superiority reached
        eval(parse(text = paste(paste0("sup_reached[", i, ",", j, "]"), "<-",
                                get(paste0("prob_sup_vec", i))[j] > Bayes_Sup[[j]][i,2])))

      }

    }

    # Add decisions from all decision rules for each arm together
    superiority_reached <- c(superiority_reached, list(sup_reached))

  }


  ########## Bayesian Two-Arm Futility Criteria ##########

  if (!is.null(Bayes_Fut1)) {
    # Get the decision rule for interim/final
    if (analysis_number == 1) {
      Bayes_Fut1 <- Bayes_Fut1[[1]]
      Bayes_Fut2 <- Bayes_Fut2[[1]]
    } else if (analysis_number == 2) {
      Bayes_Fut1 <- Bayes_Fut1[[2]]
      Bayes_Fut2 <- Bayes_Fut2[[2]]
    } else if (analysis_number == 3) {
      Bayes_Fut1 <- Bayes_Fut1[[3]]
      Bayes_Fut2 <- Bayes_Fut2[[3]]
    }

    if (endpoint_number == 1) {
      Bayes_Fut <- Bayes_Fut1
    } else if (endpoint_number == 2) {
      Bayes_Fut <- Bayes_Fut2
    }


    # Firstly initiate as many vectors as will maximum be needed (which is determined by the most complex decision rule)
    for (k in 1:max(sapply(Bayes_Fut, function(x) nrow(x)))) {
      assign(paste0("prob_fut_vec", k), rep(NA, 1))
    }

    assign(paste0("fut_reached"), matrix(NA, nrow = max(sapply(Bayes_Fut, function(x) nrow(x))), ncol = 1))

    # For each of the comparisons separately, evaluate the applying decision rules
    # j=1: Combo vs Mono, j=2: Combo vs Back, j=3: Mono vs Placebo

    for (i in 1:max(sapply(Bayes_Fut, function(x) nrow(x)))) {
      for (j in 1:1) {
        # Get which type of comparison should be conducted
        c1 <- 1
        c2 <- 2

        # Evaluate decision rule, but if non-existant, give NA
        if (is.na(Bayes_Fut[[j]][i,1])) {
          eval(parse(text = paste(paste0("prob_fut_vec", i, "[", j, "]"), "<- NA")))
        } else {
          eval(parse(text = paste(paste0("prob_fut_vec", i, "[", j, "]"), "<-",
                                  post_prob_bin(
                                    n_tot[c1], n_tot[c2],
                                    resp_tot[c1], resp_tot[c2],
                                    Bayes_Fut[[j]][i,1], beta_prior, beta_prior, beta_prior, beta_prior))))
        }

        # Check whether futility reached
        eval(parse(text = paste(paste0("fut_reached[", i, ",", j, "]"), "<-",
                                get(paste0("prob_fut_vec", i))[j] < Bayes_Fut[[j]][i,2])))
      }

    }

    # Add decisions from all decision rules for each arm together
    futility_reached <- c(futility_reached, list(fut_reached))

  }


  ########## Combine Results And Return List ##########

  # Helper Function
  check_fun <- function(x, type) {
    if (length(x) == 0 | all(sapply(x, function(y) all(is.na(y))))) {
      ret <- FALSE
    } else {
      if (type != "fut") {
        ret <- all(sapply(x, function(y) all(y, na.rm = TRUE)), na.rm = TRUE)
      } else {
        ret <- any(sapply(x, function(y) any(y, na.rm = TRUE)), na.rm = TRUE)
      }
    }
    return(ret)
  }

  # Get final results, depending on whether interim or final
  # Currently all single decisions need to be TRUE in order for the total to be TRUE

  sup <- check_fun(superiority_reached, "sup")
  fut <- check_fun(futility_reached, "fut")

  eval(parse(text = paste(paste0("res_list[[which_cohort]]$Meta$sup_interim",
                                 endpoint_number, analysis_number, "list", "<- superiority_reached"))))
  eval(parse(text = paste(paste0("res_list[[which_cohort]]$Meta$fut_interim",
                                 endpoint_number, analysis_number, "list", "<- futility_reached"))))
  eval(parse(text = paste(paste0("res_list[[which_cohort]]$Meta$sup_interim",
                                 endpoint_number, analysis_number, "<- sup"))))
  eval(parse(text = paste(paste0("res_list[[which_cohort]]$Meta$fut_interim",
                                 endpoint_number, analysis_number, "<- fut"))))


  if (fut) {
    eval(parse(text = paste(paste0("res_list[[which_cohort]]$Meta$decision_hist", endpoint_number,
                                   "[", analysis_number, "] <- 'STOP_FUT'"))))
  } else if (sup) {
    eval(parse(text = paste(paste0("res_list[[which_cohort]]$Meta$decision_hist", endpoint_number,
                                   "[", analysis_number, "] <- 'GO_SUP'"))))
  } else {
    eval(parse(text = paste(paste0("res_list[[which_cohort]]$Meta$decision_hist", endpoint_number,
                                   "[", analysis_number, "] <- 'CONTINUE'"))))
  }

  return(res_list)

}

