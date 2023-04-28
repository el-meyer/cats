#' Simulates the cohort trial.
#'
#' @param n_fin                Sample size per cohort at final
#'
#' @param cohorts_start        Number of cohorts to start the platform with
#'
#' @param rr_comb1             Response rates of treatment, histology endpoint 1
#'
#' @param rr_comb2             Response rates of treatment, histology endpoint 2
#'
#' @param rr_plac1             Response rate of the SoC, histology endpoint 1
#'
#' @param rr_plac2             Response rate of the SoC, histology endpoint 2
#'
#' @param correlation          Correlation between histology endpoints
#'
#' @param random               Should the response rates of the arms be randomly drawn from rr_exp? Default is FALSE.
#'
#' @param random_type          How should the response rates be drawn randomly? Options are:
#'
#'                             "absolute": Specify absolute response rates that will be drawn with a certain probability
#'
#'                             "risk_difference": Specify absolute response rates for placebo which will be drawn randomly,
#'                             plus specify vectors for absolute treatment effects of mono therapies over placebo and
#'                             for combo over the mono therapies.
#'
#'                             "risk_ratio": Specify absolute response rates for placebo which will be drawn randomly, plus specify vectors
#'                             for relative treatment effects of mono therapies over placebo and for combo over the mono therapies.
#'
#'                             "odds_ratios": Specify response rate for placebo, specify odds-ratios for mono therapies
#'                             (via rr_back and rr_mono) and respective probabilities.
#'                             On top, specify interaction for the combination therapy via rr_comb with prob_rr_comb.
#'
#'                             Set: odds_combo = odds_plac * or_mono1 * or_mono2 * rr_comb.
#'                             If rr_comb > 1 -> synergistic, if rr_comb = 1 -> additive. If rr_comb < 1 -> antagonistic.
#'                             Default is "NULL".
#'
#' @param prob_comb1_rr        If random == TRUE, what are the probabilities with which the elements of rr_comb1 should be drawn?
#'
#' @param prob_comb2_rr        If random == TRUE, what are the probabilities with which the elements of rr_comb2 should be drawn?
#'
#' @param prob_plac1_rr        If random == TRUE, what are the probabilities with which the elements of rr_plac1 should be drawn?
#'
#' @param prob_plac2_rr        If random == TRUE, what are the probabilities with which the elements of rr_plac2 should be drawn?
#'
#' @param stage_data           Should individual stage data be passed along? Default is TRUE
#'
#' @param cohort_random        If not NULL, indicates that new arms/cohorts should be randomly started.
#'                             For every timestep, there is a cohort_random probability that a new cohort will be started.
#'
#' @param cohort_fixed         If not NULL, fixed timesteps after which a cohort will be included
#'
#' @param cohorts_max          Maximum number of cohorts that are allowed to be added throughout the trial
#'
#' @param cohorts_sim          Maximum number of cohorts that can run simultaneously
#'
#' @param cohort_offset        Minimum number of time between adding new cohorts
#'
#' @param sr_drugs_pos         Stopping rule for successful experimental arms; Default = 1
#'
#' @param sr_pats              Stopping rule for total number of patients; Default = cohorts_max * n_fin + error term based on randomization
#'
#' @param sr_first_pos         Stopping rule for first successful cohort;
#'                             if TRUE, after first cohort was found to be successful, no further cohorts will be included
#'                             but cohorts will finish evaluating, unless other stopping rules reached prior.
#'                             Default is FALSE.
#'
#' @param sharing_type         Which backbone and placebo data should be used for arm comparisons; Default is "all".
#'                             Another option is "concurrent" or "dynamic" or "cohort".
#'
#' @param safety_prob          Probability for a random stopping after every patient
#'
#' @param missing_prob         Probability for a missing value at final (independent of treatment)
#'
#' @param accrual_type         Type of patient accrual; choices are "fixed", "poisson" or "exponential"
#'
#' @param accrual_param        Parameter used for patient accrual
#'
#' @param analysis_times       Vector of information fractions needed for first interim, second interim and final
#'
#' @param hist_lag             Time until histology outcome is observed
#'
#' @param time_trend           Additive term by which response rates increase at every time step
#'
#' @param composite            Rule for deriving the composite endpoint. By default "or", otherwise "and"
#'
#' @param ...                  Further arguments to be passed to decision function, such as decision making criteria
#'
#' @return List containing: Responses and patients on experimental and control arm, total treatment successes and failures and final p-value
#'
#' @examples
#'
#' random <- TRUE
#'
#' rr_comb1 <- 0.10
#' prob_comb1_rr <- 1
#' rr_comb2 <- 0.45
#' prob_comb2_rr <- 1
#' rr_plac1 <- 0.10
#' prob_plac1_rr <- 1
#' rr_plac2 <- 0.20
#' prob_plac2_rr <- 1
#'
#' correlation <- 0.8
#'
#' cohorts_start <- 2
#' cohorts_max <- 5
#' safety_prob <- 0
#' sharing_type <- "concurrent"
#' sr_drugs_pos <- 5
#' sr_first_pos <- FALSE
#' n_fin <- 100
#' stage_data <- TRUE
#' cohort_random <- 0.01
#' cohort_offset <- 0
#' cohorts_sim <- Inf
#' random_type <- "absolute"
#' missing_prob <- 0.2
#' cohort_fixed <- 5
#' hist_lag <- 48
#' analysis_times <- c(0.5, 0.75, 1)
#' accrual_type <- "fixed"
#' accrual_param <- 9
#' time_trend <- 0.001
#' composite <- "or"
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
#' # Comparison IA1
#' Bayes_Fut11 <- matrix(nrow = 1, ncol = 2)
#' Bayes_Fut11[1,] <- c(0.00, 0.20)
#' # Comparison IA2
#' Bayes_Fut12 <- matrix(nrow = 1, ncol = 2)
#' Bayes_Fut12[1,] <- c(0.00, 0.30)
#' # Comparison IA3
#' Bayes_Fut13 <- matrix(nrow = 1, ncol = 2)
#' Bayes_Fut13[1,] <- c(NA, NA)
#' # Endpoint 1+2
#' Bayes_Fut1 <- Bayes_Fut2 <- list(list(Bayes_Fut11), list(Bayes_Fut12), list(Bayes_Fut13))
#'
#' simulate_trial(
#'  n_fin = n_fin, random_type = random_type, composite = composite,
#'  rr_comb1 = rr_comb1, rr_comb2 = rr_comb2, rr_plac1 = rr_plac1, rr_plac2 = rr_plac2,
#'  random = random, prob_comb1_rr = prob_comb1_rr, prob_comb2_rr = prob_comb2_rr,
#'  prob_plac1_rr = prob_plac1_rr, prob_plac2_rr = prob_plac2_rr, correlation = correlation,
#'  stage_data = stage_data, cohort_random = cohort_random, cohorts_max = cohorts_max,
#'  sr_drugs_pos = sr_drugs_pos, sharing_type = sharing_type, Bayes_Fut1 = Bayes_Fut1,
#'  safety_prob = safety_prob, Bayes_Sup1 = Bayes_Sup1, Bayes_Sup2 = Bayes_Sup2,
#'  cohort_offset = cohort_offset, sr_first_pos = sr_first_pos, Bayes_Fut2 = Bayes_Fut2,
#'  missing_prob = missing_prob, cohort_fixed = cohort_fixed, accrual_type = accrual_type,
#'  accrual_param = accrual_param, hist_lag = hist_lag, analysis_times = analysis_times,
#'  time_trend = time_trend, cohorts_start = cohorts_start, cohorts_sim = cohorts_sim
#' )
#'
#' @export
simulate_trial <- function(n_fin, cohorts_start = 1, composite = "or",
                           rr_comb1, rr_plac1, rr_comb2, rr_plac2,
                           random_type = NULL, random = FALSE, correlation,
                           prob_comb1_rr = NULL, prob_plac1_rr = NULL, prob_comb2_rr = NULL,
                           prob_plac2_rr = NULL, stage_data = FALSE,
                           cohort_random = NULL, cohorts_max = 4, sr_drugs_pos = 1,
                           sr_pats = cohorts_max * (n_fin + 3 * cohorts_max), sr_first_pos = FALSE,
                           cohort_offset = 0, sharing_type = "all", safety_prob = 0, cohorts_sim = Inf,
                           missing_prob = 0, cohort_fixed = NULL, accrual_type = "fixed", accrual_param = 9,
                           hist_lag = 48, analysis_times = c(0.5, 0.75, 1), time_trend = time_trend, ...) {

  ##### Helper Functions #####

  # ------------------ Helper functions
  sample.vec <- function(x, ...) x[sample(length(x), ...)]


  # helper function check which cohort is left
  coh_left_check <- function(x) {
    if (x$Meta$decision[1] %in% c("none", "CONTINUE")
        & x$Meta$decision[2] %in% c("none", "CONTINUE")
        & x$Meta$decision[3] == "none") {
      ret <- TRUE
    } else {
      ret <- FALSE
    }
    return(ret)
  }

  # helper function check which cohort still needs to enrol patients
  coh_left_enrol_check <- function(x) {
    if (x$Meta$decision[1] %in% c("none", "CONTINUE")
        & x$Meta$decision[2] %in% c("none", "CONTINUE")
        & x$Meta$decision[3] == "none") {

      if (x$Meta$pat_enrolled < n_fin) {
        ret <- TRUE
      } else {
        ret <- FALSE
      }
    } else {
      ret <- FALSE
    }
    return(ret)
  }


  # helper function to create initial cohort
  create_cohort_initial <- function(cohorts_start, rr_comb1_vec, rr_comb2_vec, rr_plac1_vec, rr_plac2_vec) {

    res_list <-
      rep(
        list(
          list(
            Meta = list(
              decision = rep("none", 3),
              decision_hist1 = rep("none", 3),
              decision_hist2 = rep("none", 3),
              start_n = 0,
              start_time = 0,
              pat_enrolled = 0
            ),
            Arms = rep(
              list(
                list(
                  rr = NULL,
                  hist_observed = 0
                )
              ),
              2
            )
          )
        ),
        cohorts_start
      )

    arm_names <- c("Comb", "Plac")

    for (i in 1:cohorts_start) {
      names(res_list)[i] <- paste0("Cohort", i)
      names(res_list[[i]]$Arms) <- arm_names

      res_list[[i]]$Arms$Comb$rr <- matrix(c(rr_comb1_vec[i], rr_comb2_vec[i]), ncol = 2)
      res_list[[i]]$Arms$Plac$rr <- matrix(c(rr_plac1_vec[i], rr_plac2_vec[i]), ncol = 2)
    }

    return(res_list)

  }

  # helper function to create new cohort
  create_cohort_new <- function(res_list, plat_time, rr_comb1_vec, rr_comb2_vec, rr_plac1_vec, rr_plac2_vec) {

    new_list <-
        list(
          list(
            Meta = list(
              decision = rep("none", 3),
              decision_hist1 = rep("none", 3),
              decision_hist2 = rep("none", 3),
              start_n = sum(sapply(res_list, function(x) x$Meta$pat_enrolled), na.rm = T),
              start_time = plat_time,
              pat_enrolled = 0
            ),
            Arms = rep(
              list(
                list(
                  rr = 0,
                  hist_observed = 0
                )
              ),
              2
            )
          )
        )

    arm_names <- c("Comb", "Plac")

      names(new_list)[[1]] <- paste0("Cohort", length(res_list) + 1)
      names(new_list[[1]]$Arms) <- arm_names

      new_list[[1]]$Arms$Comb$rr <-
        cbind(
          pmin(
            seq(
              from = rr_comb1_vec[length(res_list) + 1],
              by = time_trend,
              length.out = plat_time)
            ,
            1
          ),
          pmin(
            seq(
              from = rr_comb2_vec[length(res_list) + 1],
              by = time_trend,
              length.out = plat_time)
            ,
            1
          )
        )

      new_list[[1]]$Arms$Plac$rr <-
        cbind(
          pmin(
            seq(
              from = rr_plac1_vec[length(res_list) + 1],
              by = time_trend,
              length.out = plat_time)
            ,
            1
          ),
          pmin(
            seq(
              from = rr_plac2_vec[length(res_list) + 1],
              by = time_trend,
              length.out = plat_time)
            ,
            1
          )
        )

    res_list <- c(res_list, new_list)

    return(res_list)

  }


  # helper function to create randomization list
  # For now, balanced to all available arms
  create_rand_list <- function (res_list) {
    # check which cohorts are left
    cohorts_left <- which(sapply(res_list, function(x) coh_left_enrol_check(x)))

    # if no cohort left, return NA, otherwise create list
    if (length(cohorts_left) == 0) {

      rand_list <- NULL

    } else {

      # get names of arms
      arm_names <- c("Comb", "Plac")

      # create randomization list
      # first column refers to cohort, second to arm
      rand_list <- matrix(nrow = length(cohorts_left) * 2, ncol = 3)

      rand_list[, 1] <- sample(rep(cohorts_left, times = 2))

      # need to go and for each cohort do a random assignment
      for (i in cohorts_left) {
        rand_list[, 2][rand_list[, 1] == i] <- sample(arm_names)
      }

      rand_list[, 3] <- 0
    }

    return(rand_list)
  }



  # helper function to check whether platform stopping rules are reached
  is_sr_reached <- function(res_list, sr_drugs_pos, sr_pats, expected) {
    ret <- 0
    positives <- sum(substring(sapply(res_list, function(x) x$Meta$decision[3]), 1, 2) == "GO")
    if (positives >= sr_drugs_pos) {
      ret <- 1
    }
    if (sr_pats < expected) {
      if (sum(sapply(res_list, function(x) total_n(x)), na.rm = T) > sr_pats) {
        ret <- 1
      }
    }
    return(ret)
  }

  # helper function to update response rates of arms
  # For now does not change
  update_rr <- function(res_list) {
    for (i in 1:length(res_list)) {
      for (j in names(res_list[[i]]$Arms)) {
        # Just add time trend
        # Make sure not larger than 1
        res_list[[i]]$Arms[[j]]$rr <-
          rbind(
            res_list[[i]]$Arms[[j]]$rr,
            pmin(utils::tail(res_list[[i]]$Arms[[j]]$rr, 1) + time_trend, 1)
          )
      }
    }
    return(res_list)
  }


  # Helper function to check whether interim milestone 1 was reached
  check_int1_milestone <- function(y, time) {

    n <- sum(sapply(y$Arms, function(x) x$hist_observed))/n_fin >= time
    # check if no first interim has been conducted yet
    new <- y$Meta$decision[1] == "none"

    return(all(n, new))

  }

  # Helper function to check whether interim milestone 2 was reached
  check_int2_milestone <- function(y, time) {

    n <- sum(sapply(y$Arms, function(x) x$hist_observed))/n_fin >= time
    # check whether first interim was conducted, but second was not
    new <- (y$Meta$decision[2] == "none") & (y$Meta$decision[1] != "none")

    return(all(n, new))

  }

  # Helper function to check whether final milestone was reached
  check_fin_milestone <- function(y, time) {

    n <- sum(sapply(y$Arms, function(x) x$hist_observed))/n_fin >= time
    # check whether first and second interim were conducted, but final was not
    new <- (y$Meta$decision[3] == "none") & (y$Meta$decision[1] != "none") & (y$Meta$decision[2] != "none")

    return(all(n, new))

  }


  # Helper function to observe outcomes
  observe_outcomes <- function(res_list, CurrentTime, hist_lag) {

    for (k in 1:length(res_list)) {
      for (l in 1:2) {

        ind_hist_obs <-
          which(
            df$Cohort == k &
              df$Arm == names(res_list[[k]]$Arms)[l] &
              df$ArrivalTime < CurrentTime - hist_lag
          )


        res_list[[k]]$Arms[[l]]$hist_observed <-
          length(df$RespHist1[ind_hist_obs])
      }
    }

    return(res_list)
  }

  # Helper function to create multinomial distribution
  fun_multnom <- function(rr_short, rr_long, correlation) {

  prob11 <-
    as.numeric(
      mvtnorm::pmvnorm(
        upper = c(stats::qnorm(rr_short), stats::qnorm(rr_long)),
        corr = matrix(
          c(1, correlation,
            correlation, 1),
          nrow = 2, ncol = 2
        )
      )
    )

  prob01 <-
    as.numeric(
      mvtnorm::pmvnorm(
        lower = c(stats::qnorm(rr_short), -Inf),
        upper = c(Inf, stats::qnorm(rr_long)),
        corr = matrix(
          c(1, correlation,
            correlation, 1),
          nrow = 2, ncol = 2
        )
      )
    )

  prob10 <-
    as.numeric(
      mvtnorm::pmvnorm(
        lower = c(-Inf, stats::qnorm(rr_long)),
        upper = c(stats::qnorm(rr_short), Inf),
        corr = matrix(
          c(1, correlation,
            correlation, 1),
          nrow = 2, ncol = 2
        )
      )
    )


  prob00 <-
    as.numeric(mvtnorm::pmvnorm(
      lower = c(stats::qnorm(rr_short), stats::qnorm(rr_long)),
      upper = c(Inf, Inf),
      corr = matrix(
        c(1, correlation,
          correlation, 1),
        nrow = 2, ncol = 2
      )
    )
    )

  return(
    c(
      p00 = prob00,
      p10 = prob10,
      p01 = prob01,
      p11 = prob11
    )
  )

}

  ###### Initialization ######

  # Get response rates for the arms
  # Check whether random experimental response rates.
  # If so, simulate, if not, then length must equal cohorts_max (including first cohort)

  if (random) {
    if (random_type == "absolute") {
      # Sample response rates for all possible arms
      rr_comb1_vec <- sample.vec(rr_comb1, cohorts_max, prob = prob_comb1_rr, replace = TRUE)
      rr_comb2_vec <- sample.vec(rr_comb2, cohorts_max, prob = prob_comb2_rr, replace = TRUE)
      rr_plac1_vec <- sample.vec(rr_plac1, cohorts_max, prob = prob_plac1_rr, replace = TRUE)
      rr_plac2_vec <- sample.vec(rr_plac2, cohorts_max, prob = prob_plac2_rr, replace = TRUE)
    }

    if (random_type == "risk_difference") {
      rr_plac1_vec <- sample.vec(rr_plac1, cohorts_max, prob = prob_plac1_rr, replace = TRUE)
      rr_plac2_vec <- sample.vec(rr_plac2, cohorts_max, prob = prob_plac2_rr, replace = TRUE)
      comb1_add    <- sample.vec(rr_comb1, cohorts_max, prob = prob_comb1_rr, replace = TRUE)
      comb2_add    <- sample.vec(rr_comb2, cohorts_max, prob = prob_comb2_rr, replace = TRUE)
      rr_comb1_vec <- pmin(rr_plac1_vec + comb1_add, 1)
      rr_comb2_vec <- pmin(rr_plac2_vec + comb2_add, 1)
    }

    if (random_type == "risk_ratio") {
      rr_plac1_vec <- sample.vec(rr_plac1, cohorts_max, prob = prob_plac1_rr, replace = TRUE)
      rr_plac2_vec <- sample.vec(rr_plac2, cohorts_max, prob = prob_plac2_rr, replace = TRUE)
      comb1_add    <- sample.vec(rr_comb1, cohorts_max, prob = prob_comb1_rr, replace = TRUE)
      comb2_add    <- sample.vec(rr_comb2, cohorts_max, prob = prob_comb2_rr, replace = TRUE)
      rr_comb1_vec <- pmin(rr_plac1_vec * comb1_add, 1)
      rr_comb2_vec <- pmin(rr_plac2_vec * comb2_add, 1)
    }

    if (random_type == "odds_ratios") {
      odds_to_rr <- function(x) {x/(1+x)}
      rr_to_odds <- function(x) {x/(1-x)}
      # get placebo response rate
      rr_plac1_vec <- sample.vec(rr_plac1, cohorts_max, prob = prob_plac1_rr, replace = TRUE)
      rr_plac2_vec <- sample.vec(rr_plac2, cohorts_max, prob = prob_plac2_rr, replace = TRUE)
      # get back, mono, comb odds ratios
      comb1_add_or <- sample.vec(rr_comb1, cohorts_max, prob = prob_comb1_rr, replace = TRUE)
      comb2_add_or <- sample.vec(rr_comb2, cohorts_max, prob = prob_comb2_rr, replace = TRUE)
      # compute mono and back odds
      odds_plac1_vec <- rr_to_odds(rr_plac1_vec)
      odds_plac2_vec <- rr_to_odds(rr_plac2_vec)
      odds_comb1_vec <- odds_plac1_vec * comb1_add_or
      odds_comb2_vec <- odds_plac2_vec * comb2_add_or
      # transfer odds to rr
      rr_comb1_vec <- odds_to_rr(odds_comb1_vec)
      rr_comb2_vec <- odds_to_rr(odds_comb2_vec)
    }

  } else {
    rr_comb1_vec <- rep(rr_comb1, cohorts_max)
    rr_comb2_vec <- rep(rr_comb2, cohorts_max)
    rr_plac1_vec <- rep(rr_plac1, cohorts_max)
    rr_plac2_vec <- rep(rr_plac2, cohorts_max)
  }


  # dummy to indicate trial stop
  trial_stop <- 0

  # dummy for timestamp for first success
  first_success <- -1

  # Variable measuring time since last cohort was added
  last_cohort_time <- 0

  # Initialize res_list
  res_list <- create_cohort_initial(cohorts_start, rr_comb1_vec, rr_comb2_vec, rr_plac1_vec, rr_plac2_vec)

  Total_N_Vector <- NULL

  # initialize platform time
  plat_time <- 0

  # Initialize new Patient vector and vector of exact arrival times
  new_pats <- NULL
  pats_arrival_times <- NULL

  # Initialize empty data frame
  cols <- c("PatID", "ArrivalTime", "Cohort", "Arm", "RespHist1", "RespHist2", "HistMissing")
  df <- matrix(nrow = 0, ncol = length(cols))
  colnames(df) <- cols
  df <- as.data.frame(df)

  # create initial randomization list
  rand_list <- create_rand_list(res_list)

  ##### Running Simulations #####
  while (!trial_stop) {

    plat_time <- plat_time + 1

    ##### Accrue new patients ######

    if (accrual_type == "fixed") {
      new_pats <- c(new_pats, accrual_param)
    }

    if (accrual_type == "poisson") {
      new_pats <- c(new_pats, stats::rpois(1, accrual_param))
    }

    if (accrual_type == "exponential") {
      new_pats <- c(new_pats, round(accrual_param ^ (plat_time)))
    }

    new_arrival_times <- sort(stats::runif(new_pats[plat_time])) + plat_time - 1
    pats_arrival_times <- c(pats_arrival_times, new_arrival_times)

    # Loop over patient arrival times

    for (i in new_arrival_times) {

      # only add patients to cohorts, if any arm is still enrolling, i.e. randomization list is not NA
      if (!is.null(rand_list)) {

        # get index of next line in randomization list and if there are only used entries, create a new one
        if (suppressWarnings((min(which(rand_list[, 3] == 0)))) == Inf) {
          rand_list <- create_rand_list(res_list)
          index <- 1
          rand_list[1, 3] <- 1
        } else {
          index <- min(which(rand_list[, 3] == 0))
          rand_list[index, 3] <- 1
        }

        # Cohort can be found in the first column of randomization list
        new_pat_cohort <- as.numeric(rand_list[index, 1])

        # Name of Arm can be found in the second column of the randomization list
        new_pat_arm <- rand_list[index, 2]

        # add one more patient to the cohort
        res_list[[new_pat_cohort]]$Meta$pat_enrolled <- res_list[[new_pat_cohort]]$Meta$pat_enrolled + 1

        # if final sample size is reached, create new randomization list
        if (res_list[[new_pat_cohort]]$Meta$pat_enrolled == n_fin) {
          rand_list <- create_rand_list(res_list)
        }

        # Simulate multinomial outcome of patient
        rr_hist1 <- res_list[[new_pat_cohort]]$Arms[[new_pat_arm]]$rr[plat_time, 1]
        rr_hist2 <- res_list[[new_pat_cohort]]$Arms[[new_pat_arm]]$rr[plat_time, 2]
        new_probs <- fun_multnom(rr_hist1, rr_hist2, correlation = correlation)
        draw <- t(stats::rmultinom(1, 1, new_probs))
        resp_hist1 <- 0
        resp_hist2 <- 0
        if (draw[1,2] == 1) {
          resp_hist1 <- resp_hist1 + 1
        }
        if (draw[1,3] == 1) {
          resp_hist2 <- resp_hist2 + 1
        }
        if (draw[1,4] == 1) {
          resp_hist1 <- resp_hist1 + 1
          resp_hist2 <- resp_hist2 + 1
        }


        # create data of new patient
        single_pat <-
          data.frame(
            PatID       = nrow(df) + 1, # Patient ID is one more than what was there previously
            ArrivalTime = i, # current exact time
            Cohort      = new_pat_cohort,
            Arm         = new_pat_arm,
            RespHist1   = resp_hist1,
            RespHist2   = resp_hist2,
            HistMissing = sample(0:1, 1, prob = c(1 - missing_prob, missing_prob)) # sample whether hist endpoint will be missing
          )

        # Add patient to existing data frame
        df <- rbind(df, single_pat)
      }

      # Update the Meta information for every arm/cohort (i.e. how many outcomes are observed)
      res_list <- observe_outcomes(res_list, i, hist_lag)

      # Check whether any analyses should be conducted now

      ind_interim1 <- which(sapply(res_list, function(x) check_int1_milestone(x, time = analysis_times[1])))
      ind_interim2 <- which(sapply(res_list, function(x) check_int2_milestone(x, time = analysis_times[2])))
      ind_final   <- which(sapply(res_list, function(x) check_fin_milestone(x, time = analysis_times[3])))

      if (length(ind_interim1) > 0) {
        ############# conduct interim analysis
        for (j in ind_interim1) {

          res_list <-
            make_decision_trial(
              res_list         = res_list,
              which_cohort     = j,
              analysis_number  = 1,
              analysis_time    = i,
              hist_lag         = hist_lag,
              dataset          = df,
              sharing_type     = sharing_type,
              endpoint_number  = 1,
              ...
            )

          res_list <-
            make_decision_trial(
              res_list         = res_list,
              which_cohort     = j,
              analysis_number  = 1,
              analysis_time    = i,
              hist_lag         = hist_lag,
              dataset          = df,
              sharing_type     = sharing_type,
              endpoint_number  = 2,
              ...
            )

          # Get combined decision
          if (composite == "or") {

            # if any of the two endpoints is successful, declare success
            if (res_list[[j]]$Meta$decision_hist1[1] == "GO_SUP"
                | res_list[[j]]$Meta$decision_hist2[1] == "GO_SUP") {
              res_list[[j]]$Meta$decision[1] <- "GO_SUP"
            } else
            # if both of the two endpoints are futile, declare futility
            if (res_list[[j]]$Meta$decision_hist1[1] == "STOP_FUT"
                & res_list[[j]]$Meta$decision_hist2[1] == "STOP_FUT") {
              res_list[[j]]$Meta$decision[1] <- "STOP_FUT"
            } else {
              res_list[[j]]$Meta$decision[1] <- "CONTINUE"
            }
          } else if (composite == "and") {
          # opposite of before
            if (res_list[[j]]$Meta$decision_hist1[1] == "GO_SUP"
                & res_list[[j]]$Meta$decision_hist2[1] == "GO_SUP") {
              res_list[[j]]$Meta$decision[1] <- "GO_SUP"
            } else
              if (res_list[[j]]$Meta$decision_hist1[1] == "STOP_FUT"
                  | res_list[[j]]$Meta$decision_hist2[1] == "STOP_FUT") {
                res_list[[j]]$Meta$decision[1] <- "STOP_FUT"
              } else {
                res_list[[j]]$Meta$decision[1] <- "CONTINUE"
              }
          }

          # What happens at successful interim
          if (res_list[[j]]$Meta$decision[1] == "GO_SUP") {
            res_list[[j]]$Meta$decision[2] <- "GO_SUP"
            res_list[[j]]$Meta$decision[3] <- "GO_SUP"
            if (first_success == -1) {
              first_success <- plat_time
            }
            rand_list <- create_rand_list(res_list)
          }

          # What happens at unsuccessful interim
          if (res_list[[j]]$Meta$decision[1] == "STOP_FUT") {
            res_list[[j]]$Meta$decision[2] <- "STOP_FUT"
            res_list[[j]]$Meta$decision[3] <- "STOP_FUT"
            rand_list <- create_rand_list(res_list)
          }

        }
      }

      if (length(ind_interim2) > 0) {
        ############# conduct interim analysis
        for (j in ind_interim2) {

          res_list <-
            make_decision_trial(
              res_list         = res_list,
              which_cohort     = j,
              analysis_number  = 2,
              analysis_time    = i,
              hist_lag         = hist_lag,
              dataset          = df,
              sharing_type     = sharing_type,
              endpoint_number  = 1,
              ...
            )

          res_list <-
            make_decision_trial(
              res_list         = res_list,
              which_cohort     = j,
              analysis_number  = 2,
              analysis_time    = i,
              hist_lag         = hist_lag,
              dataset          = df,
              sharing_type     = sharing_type,
              endpoint_number  = 2,
              ...
            )

          # Get combined decision
          if (composite == "or") {

            # if any of the two endpoints is successful, declare success
            if (res_list[[j]]$Meta$decision_hist1[2] == "GO_SUP"
                | res_list[[j]]$Meta$decision_hist2[2] == "GO_SUP") {
              res_list[[j]]$Meta$decision[2] <- "GO_SUP"
            } else
              # if both of the two endpoints are futile, declare futility
              if (res_list[[j]]$Meta$decision_hist1[2] == "STOP_FUT"
                  & res_list[[j]]$Meta$decision_hist2[2] == "STOP_FUT") {
                res_list[[j]]$Meta$decision[2] <- "STOP_FUT"
              } else {
                res_list[[j]]$Meta$decision[2] <- "CONTINUE"
              }
          } else if (composite == "and") {
            # opposite of before
            if (res_list[[j]]$Meta$decision_hist1[2] == "GO_SUP"
                & res_list[[j]]$Meta$decision_hist2[2] == "GO_SUP") {
              res_list[[j]]$Meta$decision[2] <- "GO_SUP"
            } else
              if (res_list[[j]]$Meta$decision_hist1[2] == "STOP_FUT"
                  | res_list[[j]]$Meta$decision_hist2[2] == "STOP_FUT") {
                res_list[[j]]$Meta$decision[2] <- "STOP_FUT"
              } else {
                res_list[[j]]$Meta$decision[2] <- "CONTINUE"
              }
          }

          # What happens at successful interim
          if (res_list[[j]]$Meta$decision[2] == "GO_SUP") {
            res_list[[j]]$Meta$decision[3] <- "GO_SUP"
            if (first_success == -1) {
              first_success <- plat_time
            }
            rand_list <- create_rand_list(res_list)
          }

          # What happens at unsuccessful interim
          if (res_list[[j]]$Meta$decision[2] == "STOP_FUT") {
            res_list[[j]]$Meta$decision[3] <- "STOP_FUT"
            rand_list <- create_rand_list(res_list)
          }

        }
      }

      if (length(ind_final) > 0) {
        ############# conduct interim analysis
        for (j in ind_final) {

          res_list <-
            make_decision_trial(
              res_list         = res_list,
              which_cohort     = j,
              analysis_number  = 3,
              analysis_time    = i,
              hist_lag         = hist_lag,
              dataset          = df,
              sharing_type     = sharing_type,
              endpoint_number  = 1,
              ...
            )

          res_list <-
            make_decision_trial(
              res_list         = res_list,
              which_cohort     = j,
              analysis_number  = 3,
              analysis_time    = i,
              hist_lag         = hist_lag,
              dataset          = df,
              sharing_type     = sharing_type,
              endpoint_number  = 2,
              ...
            )

          # Get combined decision
          if (composite == "or") {

            # if any of the two endpoints is successful, declare success
            if (res_list[[j]]$Meta$decision_hist1[3] == "GO_SUP"
                | res_list[[j]]$Meta$decision_hist2[3] == "GO_SUP") {
              res_list[[j]]$Meta$decision[3] <- "GO_SUP"
            } else
              # if both of the two endpoints are futile, declare futility
              if (res_list[[j]]$Meta$decision_hist1[3] == "STOP_FUT"
                  & res_list[[j]]$Meta$decision_hist2[3] == "STOP_FUT") {
                res_list[[j]]$Meta$decision[3] <- "STOP_FUT"
              } else {
                res_list[[j]]$Meta$decision[3] <- "STOP_N"
              }
          } else if (composite == "and") {
            # opposite of before
            if (res_list[[j]]$Meta$decision_hist1[3] == "GO_SUP"
                & res_list[[j]]$Meta$decision_hist2[3] == "GO_SUP") {
              res_list[[j]]$Meta$decision[3] <- "GO_SUP"
            } else
              if (res_list[[j]]$Meta$decision_hist1[3] == "STOP_FUT"
                  | res_list[[j]]$Meta$decision_hist2[3] == "STOP_FUT") {
                res_list[[j]]$Meta$decision[3] <- "STOP_FUT"
              } else {
                res_list[[j]]$Meta$decision[3] <- "STOP_N"
              }
          }

          if (res_list[[j]]$Meta$decision[3] == "GO_SUP") {
            if (first_success == -1) {
              first_success <- plat_time
            }
            rand_list <- create_rand_list(res_list)
          }
        }
      }

      # random safety stopping
      random_stop <- sample(0:1, 1, prob = c(1 - safety_prob, safety_prob))
      if (random_stop) {
        if (res_list[[new_patcohort]]$Meta$decision[1] == "none")  {res_list[[new_patcohort]]$Meta$decision[1] <- "STOP_SAFETY"}
        if (res_list[[new_patcohort]]$Meta$decision[2] == "none")  {res_list[[new_patcohort]]$Meta$decision[2] <- "STOP_SAFETY"}
        res_list[[new_patcohort]]$Meta$decision[3] <- "STOP_SAFETY"
        rand_list <- create_rand_list(res_list)
      }

      # If all cohorts are stopped, stop trial, but only if maximum number of cohorts are not reached
      # Otherwise break loop

      if (!any(sapply(res_list, function(x) x$Meta$decision[3]) %in% c("none", "CONTINUE")) &
          length(res_list) == cohorts_max) {
        trial_stop <- 1
      }

      if (!any(sapply(res_list, function(x) x$Meta$decision[3]) %in% c("none", "CONTINUE"))) {
        break
      }
    }

    # see whether trial has stopped
    if (trial_stop == 1) {
      break
    }

    # If trial continues, check whether new cohorts should be included
    coh_add <- 0

    # Check whether any new cohorts should be added to the trial
    # is maximum number of cohorts reached
    if (length(res_list) < cohorts_max) {

      # Check whether no cohort is active, then one is added immediately
      if (!any(sapply(res_list, function(x) x$Meta$decision[3]) %in% c("none", "CONTINUE"))) {
        coh_add <- 1
        } else {

        # is maximum number of parallel cohorts reached
        if (sum(sapply(res_list, function(x) coh_left_check(x))) < cohorts_sim) {
          # is cohort offset satisfied
          if (plat_time - last_cohort_time > cohort_offset) {
            # check random adding
            if (!is.null(cohort_random)) {
              coh_add <- coh_add + sample(0:1, 1, prob = c(1 - cohort_random, cohort_random))
            }
            # check fixed schedule adding
            if (!is.null(cohort_fixed)) {
              if (plat_time - last_cohort_time >= cohort_fixed) {
                coh_add <- coh_add +1
              }
            }
          }
        }

      }
    }

    # make sure that number of arms to add is bounded above by cohorts_max and cohorts_sim
    coh_add <- min(coh_add, cohorts_max - length(res_list), cohorts_sim - sum(sapply(res_list, function(x) coh_left_check(x))))

    # Add new cohorts and create new randomization list immediately
    if (coh_add > 0) {
      for (c in 1:coh_add) {
        res_list <-
          create_cohort_new(
            res_list,
            plat_time,
            rr_comb1_vec,
            rr_comb2_vec,
            rr_plac1_vec,
            rr_plac2_vec
          )
      }

      last_cohort_time <- plat_time
      rand_list <- create_rand_list(res_list)
    }

    # if trial continues, update response rates
    res_list <- update_rr(res_list)

  }

  # Define truth via:
  # Is RR1 > RR2

  truth <- rep(NA, length(res_list))


  for (i in 1:length(res_list)) {
    if (composite == "or") {
      truth[i] <-
        (res_list[[i]]$Arms$Comb$rr[1,1] > res_list[[i]]$Arms$Plac$rr[1,1])|
        (res_list[[i]]$Arms$Comb$rr[1,2] > res_list[[i]]$Arms$Plac$rr[1,2])
    } else if (composite == "and") {
      truth[i] <-
        (res_list[[i]]$Arms$Comb$rr[1,1] > res_list[[i]]$Arms$Plac$rr[1,1])&
        (res_list[[i]]$Arms$Comb$rr[1,2] > res_list[[i]]$Arms$Plac$rr[1,2])
    }
  }

  # Get final experimental response rates over time
  rr_comb1_final <- sapply(res_list, function(x) x$Arms$Comb$rr[,1])
  rr_comb2_final <- sapply(res_list, function(x) x$Arms$Comb$rr[,2])
  rr_plac1_final <- sapply(res_list, function(x) x$Arms$Plac$rr[,1])
  rr_plac2_final <- sapply(res_list, function(x) x$Arms$Plac$rr[,2])

  # Time and patients enrolled previous to first success
  # Get time stamp of first success and then compute numbers

  if (first_success > 0) {

    time_to_first_success <- first_success
    pat_to_first_success <- nrow(df[which(df$ArrivalTime < first_success),])

  } else {

    time_to_first_success <- NA
    pat_to_first_success <- NA

  }

  # Check which decisions were correct positives, false positives, etc.
  cp <- sum(substring(sapply(res_list, function(x) x$Meta$decision[3]), 1, 2) == "GO" &  truth)
  fp <- sum(substring(sapply(res_list, function(x) x$Meta$decision[3]), 1, 2) == "GO" & !truth)
  cn <- sum(substring(sapply(res_list, function(x) x$Meta$decision[3]), 1, 2) == "ST" & !truth)
  fn <- sum(substring(sapply(res_list, function(x) x$Meta$decision[3]), 1, 2) == "ST" &  truth)

  # Prepare return list
  ret <- list(
    Decision               = sapply(res_list, function(x) x$Meta$decision),
    Start_N                = sapply(res_list, function(x) x$Meta$start_n),
    Start_Time             = sapply(res_list, function(x) x$Meta$start_time),
    RR_Comb1               = rr_comb1_final,
    RR_Comb2               = rr_comb2_final,
    RR_Plac1               = rr_plac1_final,
    RR_Plac2               = rr_plac2_final,
    N_Cohorts              = length(res_list),
    Final_N_Cohort         = sapply(res_list, function(x) x$Meta$pat_enrolled),
    Final_N_Cohort_Trial   = mean(sapply(res_list, function(x) x$Meta$pat_enrolled)),
    Total_N                = sum(sapply(res_list, function(x) x$Meta$pat_enrolled)),
    Total_Time             = plat_time,
    Time_First_Suc         = time_to_first_success,
    Pat_First_Suc          = pat_to_first_success,
    Pat_Arrival_Times      = pats_arrival_times,
    Unenrolled_Pats        = length(pats_arrival_times) - sum(sapply(res_list, function(x) x$Meta$pat_enrolled)),
    TP                     = cp,
    FP                     = fp,
    TN                     = cn,
    FN                     = fn,
    FDR_Trial              = ifelse(!is.na(fp/(cp + fp)), fp/(cp + fp), NA),
    PTP_Trial              = ifelse(!is.na(cp/(cp + fn)), cp/(cp + fn), NA),
    PTT1ER_Trial           = ifelse(!is.na(fp/(fp + cn)), fp/(fp + cn), NA),
    any_P                  = as.numeric((cp + fp) > 0),

    Intx1_GO                = sum(sapply(res_list, function(x) x$Meta$decision[1] == "GO_SUP"), na.rm = TRUE),
    Intx1_STOP              = sum(sapply(res_list, function(x) x$Meta$decision[1] == "STOP_FUT"), na.rm = TRUE),

    Intx2_GO                = sum(sapply(res_list, function(x) x$Meta$decision[2] == "GO_SUP"), na.rm = TRUE),
    Intx2_STOP              = sum(sapply(res_list, function(x) x$Meta$decision[2] == "STOP_FUT"), na.rm = TRUE),

    Intx3_GO                = sum(sapply(res_list, function(x) x$Meta$decision[3] == "GO_SUP"), na.rm = TRUE),
    Intx3_STOP              = sum(sapply(res_list, function(x) x$Meta$decision[3] == "STOP_FUT"), na.rm = TRUE),

    # Int11_GO                = sum(sapply(res_list, function(x) x$Meta$sup_interim11), na.rm = TRUE),
    # Int11_STOP              = sum(sapply(res_list, function(x) x$Meta$fut_interim11), na.rm = TRUE),
    # Int21_GO                = sum(sapply(res_list, function(x) x$Meta$sup_interim21), na.rm = TRUE),
    # Int21_STOP              = sum(sapply(res_list, function(x) x$Meta$fut_interim21), na.rm = TRUE),

    # Int12_GO                = sum(sapply(res_list, function(x) x$Meta$sup_interim12), na.rm = TRUE),
    # Int12_STOP              = sum(sapply(res_list, function(x) x$Meta$fut_interim12), na.rm = TRUE),
    # Int22_GO                = sum(sapply(res_list, function(x) x$Meta$sup_interim22), na.rm = TRUE),
    # Int22_STOP              = sum(sapply(res_list, function(x) x$Meta$fut_interim22), na.rm = TRUE),

    Safety_STOP            = sum(sapply(res_list, function(x) (x$Meta$decision[3] == "STOP_SAFETY")), na.rm = TRUE)
  )

  if (stage_data) {
    ret <- list(Trial_Overview = ret, Stage_Data = res_list, Pat_Data = df)
  } else {
    ret <- list(Trial_Overview = ret)
  }

  return(ret)

}
