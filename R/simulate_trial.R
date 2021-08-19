#' Simulates the cohort trial.
#'
#' @param n_int                Sample size per cohort to conduct interim analysis
#'
#' @param n_fin                Sample size per cohort at final
#'
#' @param cohorts_start        Number of cohorts to start the platform with
#'
#' @param design_type          Either "combination" or "doses"
#'
#' @param arms_per_cohort      An integer between 2 and 4.
#'                             If 2, then only comb and SoC.
#'                             If 3, then only Comb, Mono and SoC.
#'                             If 4, then all four arms.
#'
#' @param rr_comb              Response rates of combination therapies // high dose
#'
#' @param rr_mono              Response rate of mono therapies // medium dose
#'
#' @param rr_back              Response rates of backbone arms  // low dose
#'
#' @param rr_plac              Response rate of the placebo // SoC
#'
#' @param rr_transform         Function transforming all the above response rates to a vector of four
#'                             probabilities for the multinomial simulation
#'                             First element is probability of both failures.
#'                             Second element is probability of biomarker success and histology failure.
#'                             Third element is probability of biomarker failure and histology success.
#'                             Fourth element is probability of both success.
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
#' @param prob_comb_rr         If random == TRUE, what are the probabilities with which the elements of rr_comb should be drawn?
#'
#' @param prob_mono_rr         If random == TRUE, what are the probabilities with which the elements of rr_mono should be drawn?
#'
#' @param prob_back_rr         If random == TRUE, what are the probabilities with which the elements of rr_back should be drawn?
#'
#' @param prob_plac_rr         If random == TRUE, what are the probabilities with which the elements of rr_plac should be drawn?
#'
#' @param prob_rr_transform    If random == TRUE, what are the probabilities with which the elements of rr_transform should be drawn?
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
#' @param bio_lag              Time until biomarker outcome is observed
#'
#' @param hist_lag             Time until histology outcome is observed
#'
#' @param time_trend           Additive term by which response rates increase at every time step
#'
#' @param ...                  Further arguments to be passed to decision function, such as decision making criteria
#'
#' @return List containing: Responses and patients on experimental and control arm, total treatment successes and failures and final p-value
#'
#' @examples
#'
#' random <- TRUE
#'
#' rr_comb <- c(0.20, 0.25, 0.3)
#' prob_comb_rr <- c(0.4, 0.4, 0.2)
#' rr_plac <- c(0.10, 0.12, 0.14)
#' prob_plac_rr <- c(0.25, 0.5, 0.25)
#' rr_back <- 1
#' prob_back_rr <- 1
#' rr_mono <- 1
#' prob_mono_rr <- 1
#'
#' rr_transform <- list(function(x) {return(c(0.85*(1 - x), (1-0.85)*(1-x), (1-0.85)*x, 0.85*x))})
#' prob_rr_transform <- 1
#'
#' cohorts_start <- 2
#' design_type <- "combination"
#' arms_per_cohort <- 2
#' cohorts_max <- 5
#' safety_prob <- 0
#' sharing_type <- "concurrent"
#' sr_drugs_pos <- 5
#' sr_first_pos <- FALSE
#' n_int <- 50
#' n_fin <- 100
#' stage_data <- TRUE
#' cohort_random <- 0.01
#' cohort_offset <- 0
#' cohorts_sim <- Inf
#' random_type <- "absolute"
#' missing_prob <- 0.2
#' cohort_fixed <- 5
#' bio_lag <- 12
#' hist_lag <- 48
#' accrual_type <- "fixed"
#' accrual_param <- 9
#' time_trend <- 0.001
#'
#' # Vergleich Combo vs Plac
#' Bayes_Sup1 <- matrix(nrow = 1, ncol = 2)
#' Bayes_Sup1[1,] <- c(0.00, 0.95)
#' Bayes_Sup <- list(list(Bayes_Sup1),
#'                   list(Bayes_Sup1))
#'
#' a <- simulate_trial(
#' n_int = n_int, n_fin = n_fin, random_type = random_type,
#' rr_comb = rr_comb, rr_mono = rr_mono, rr_back = rr_back, rr_plac = rr_plac,
#' rr_transform = rr_transform, random = random, prob_comb_rr = prob_comb_rr,
#' prob_mono_rr = prob_mono_rr, prob_back_rr = prob_back_rr, prob_plac_rr = prob_plac_rr,
#' stage_data = stage_data, cohort_random = cohort_random, cohorts_max = cohorts_max,
#' sr_drugs_pos = sr_drugs_pos, sharing_type = sharing_type,design_type = design_type,
#' safety_prob = safety_prob, Bayes_Sup = Bayes_Sup, prob_rr_transform = prob_rr_transform,
#' cohort_offset = cohort_offset, sr_first_pos = sr_first_pos, arms_per_cohort = arms_per_cohort,
#' missing_prob = missing_prob, cohort_fixed = cohort_fixed, accrual_type = accrual_type,
#' accrual_param = accrual_param, bio_lag = bio_lag, hist_lag = hist_lag,
#' time_trend = time_trend, cohorts_start = cohorts_start, cohorts_sim = cohorts_sim
#' )
#'
#' @export
simulate_trial <- function(n_int, n_fin, cohorts_start = 1, design_type, arms_per_cohort,
                           rr_comb, rr_mono = 1, rr_back = 1, rr_plac,
                           rr_transform, random_type = NULL, trial_struc = "all_plac", random = FALSE,
                           prob_comb_rr = NULL, prob_mono_rr = c(1), prob_back_rr = c(1),
                           prob_plac_rr = NULL, prob_rr_transform = prob_rr_transform, stage_data = FALSE,
                           cohort_random = NULL, cohorts_max = 4, sr_drugs_pos = 1,
                           sr_pats = cohorts_max * (n_fin + 3 * cohorts_max), sr_first_pos = FALSE,
                           cohort_offset = 0, sharing_type = "all", safety_prob = 0, cohorts_sim = Inf,
                           missing_prob = 0, cohort_fixed = NULL, accrual_type = "fixed", accrual_param = 9,
                           bio_lag = 12, hist_lag = 48, time_trend = time_trend, ...) {

  ##### Helper Functions #####

  # ------------------ Helper functions
  sample.vec <- function(x, ...) x[sample(length(x), ...)]


  # helper function check which cohort is left
  coh_left_check <- function(x) {
    if (x$Meta$decision[1] %in% c("none", "PROMISING", "CONTINUE") & x$Meta$decision[2] == "none") {
      ret <- TRUE
    } else {
      ret <- FALSE
    }
    return(ret)
  }

  # helper function check which cohort still needs to enrol patients
  coh_left_enrol_check <- function(x) {
    if (x$Meta$decision[1] %in% c("none", "PROMISING", "CONTINUE") & x$Meta$decision[2] == "none") {
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
  create_cohort_initial <- function(cohorts_start, design_type, arms_per_cohort,
                                    rr_comb_vec, rr_mono_vec, rr_back_vec, rr_plac_vec) {

    res_list <-
      rep(
        list(
          list(
            Meta = list(
              decision = rep("none", 2),
              start_n = 0,
              start_time = 0,
              pat_enrolled = 0
            ),
            Arms = rep(
              list(
                list(
                  rr = NULL,
                  bio_observed = 0,
                  hist_observed = 0
                )
              ),
              arms_per_cohort
            )
          )
        ),
        cohorts_start
      )

    arm_names <- c("Comb", "Mono", "Back", "Plac")[c(1:(arms_per_cohort - 1), 4)]

    for (i in 1:cohorts_start) {
      names(res_list)[i] <- paste0("Cohort", i)
      names(res_list[[i]]$Arms) <- arm_names

      res_list[[i]]$Arms$Comb$rr <- rr_comb_vec[i]
      res_list[[i]]$Arms$Plac$rr <- rr_plac_vec[i]

      if (arms_per_cohort > 2) {
        res_list[[i]]$Arms$Mono$rr <- rr_mono_vec[i]
      }

      if (arms_per_cohort > 3) {
        res_list[[i]]$Arms$Back$rr <- rr_back_vec[i]
      }
    }

    return(res_list)

  }

  # helper function to create new cohort
  create_cohort_new <- function(res_list, plat_time, design_type, arms_per_cohort,
                                rr_comb_vec, rr_mono_vec, rr_back_vec, rr_plac_vec) {

    new_list <-
        list(
          list(
            Meta = list(
              decision = rep("none", 2),
              start_n = sum(sapply(res_list, function(x) x$Meta$pat_enrolled), na.rm = T),
              start_time = plat_time,
              pat_enrolled = 0
            ),
            Arms = rep(
              list(
                list(
                  rr = 0,
                  bio_observed = 0,
                  hist_observed = 0
                )
              ),
              arms_per_cohort
            )
          )
        )

    arm_names <- c("Comb", "Mono", "Back", "Plac")[c(1:(arms_per_cohort - 1), 4)]

      names(new_list)[[1]] <- paste0("Cohort", length(res_list) + 1)
      names(new_list[[1]]$Arms) <- arm_names

      new_list[[1]]$Arms$Comb$rr <-
        pmin(
          seq(
           from = rr_comb_vec[length(res_list) + 1],
           by = time_trend,
           length.out = plat_time)
          ,
          1
        )
       new_list[[1]]$Arms$Plac$rr <-
        pmin(
          seq(
            from = rr_plac_vec[length(res_list) + 1],
            by = time_trend,
            length.out = plat_time)
          ,
          1
         )

      if (arms_per_cohort > 2) {
        new_list[[1]]$Arms$Mono$rr <-
          pmin(
            seq(
              from = rr_mono_vec[length(res_list) + 1],
              by = time_trend,
              length.out = plat_time
            ),
            1
            )
      }

      if (arms_per_cohort > 3) {
        new_list[[1]]$Arms$Back$rr <-
          pmin(
            seq(
              from = rr_back_vec[length(res_list) + 1],
              by = time_trend,
              length.out = plat_time
              )
            ,
            1
          )
      }

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
      arm_names <- c("Comb", "Mono", "Back", "Plac")[c(1:(arms_per_cohort - 1), 4)]

      # create randomization list
      # first column refers to cohort, second to arm
      rand_list <- matrix(nrow = length(cohorts_left) * arms_per_cohort, ncol = 3)

      rand_list[, 1] <- sample(rep(cohorts_left, times = arms_per_cohort))

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
    positives <- sum(substring(sapply(res_list, function(x) x$Meta$decision[2]), 1, 2) == "GO")
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
        res_list[[i]]$Arms[[j]]$rr <- pmin(c(res_list[[i]]$Arms[[j]]$rr, tail(res_list[[i]]$Arms[[j]]$rr, 1) + time_trend), 1)
      }
    }
    return(res_list)
  }


  # Helper function to check whether interim milestone was reached
  check_int_milestone <- function(y, n_int) {

    n <- sum(sapply(y$Arms, function(x) x$bio_observed)) >= n_int
    new <- y$Meta$decision[1] == "none"

    return(all(n, new))

  }

  # Helper function to check whether final milestone was reached
  check_fin_milestone <- function(y, n_fin) {

    n <- sum(sapply(y$Arms, function(x) x$hist_observed)) >= n_fin
    new <- y$Meta$decision[2] == "none"

    return(all(n, new))

  }

  # Helper function to observe outcomes
  observe_outcomes <- function(res_list, CurrentTime, bio_lag, hist_lag) {

    for (k in 1:length(res_list)) {
      for (l in 1:arms_per_cohort) {

        ind_bio_obs <-
          which(
            df$Cohort == k &
              df$Arm == names(res_list[[k]]$Arms)[l] &
              df$ArrivalTime < CurrentTime - bio_lag
          )

        res_list[[k]]$Arms[[l]]$bio_observed <-
          length(df$RespBio[ind_bio_obs])

        ind_hist_obs <-
          which(
            df$Cohort == k &
              df$Arm == names(res_list[[k]]$Arms)[l] &
              df$ArrivalTime < CurrentTime - hist_lag
          )


        res_list[[k]]$Arms[[l]]$hist_observed <-
          length(df$RespHist[ind_hist_obs])
      }
    }

    return(res_list)
  }

  ###### Initialization ######

  # Get response rates for the arms
  # Check whether random experimental response rates.
  # If so, simulate, if not, then length must equal cohorts_max (including first cohort)

  if (random) {
    if (random_type == "absolute") {
      # Sample response rates for all possible arms
      rr_comb_vec <- sample.vec(rr_comb, cohorts_max, prob = prob_comb_rr, replace = TRUE)
      rr_back_vec <- sample.vec(rr_back, cohorts_max, prob = prob_back_rr, replace = TRUE)
      rr_mono_vec <- sample.vec(rr_mono, cohorts_max, prob = prob_mono_rr, replace = TRUE)
      rr_plac_vec <- sample.vec(rr_plac, cohorts_max, prob = prob_plac_rr, replace = TRUE)
    }

    if (random_type == "risk_difference") {
      rr_plac_vec <- sample.vec(rr_plac, cohorts_max, prob = prob_plac_rr, replace = TRUE)
      back_add    <- sample.vec(rr_back, cohorts_max, prob = prob_back_rr, replace = TRUE)
      mono_add    <- sample.vec(rr_mono, cohorts_max, prob = prob_mono_rr, replace = TRUE)
      comb_add    <- sample.vec(rr_comb, cohorts_max, prob = prob_comb_rr, replace = TRUE)
      rr_back_vec <- pmin(rr_plac_vec + back_add, 1)
      rr_mono_vec <- pmin(rr_plac_vec + mono_add, 1)
      rr_comb_vec <- pmin(rr_plac_vec + comb_add, 1)
    }

    if (random_type == "risk_ratio") {
      rr_plac_vec <- sample.vec(rr_plac, cohorts_max, prob = prob_plac_rr, replace = TRUE)
      back_add    <- sample.vec(rr_back, cohorts_max, prob = prob_back_rr, replace = TRUE)
      mono_add    <- sample.vec(rr_mono, cohorts_max, prob = prob_mono_rr, replace = TRUE)
      comb_add    <- sample.vec(rr_comb, cohorts_max, prob = prob_comb_rr, replace = TRUE)
      rr_back_vec <- pmin(rr_plac_vec * back_add, 1)
      rr_mono_vec <- pmin(rr_plac_vec * mono_add, 1)
      rr_comb_vec <- pmin(rr_plac_vec * comb_add, 1)
    }

    if (random_type == "odds_ratios") {
      odds_to_rr <- function(x) {x/(1+x)}
      rr_to_odds <- function(x) {x/(1-x)}
      # get placebo response rate
      rr_plac_vec <- sample.vec(rr_plac, cohorts_max, prob = prob_plac_rr, replace = TRUE)
      # get back, mono, comb odds ratios
      back_add_or <- sample.vec(rr_back, cohorts_max, prob = prob_back_rr, replace = TRUE)
      mono_add_or <- sample.vec(rr_mono, cohorts_max, prob = prob_mono_rr, replace = TRUE)
      comb_add_or <- sample.vec(rr_comb, cohorts_max, prob = prob_comb_rr, replace = TRUE)
      # compute mono and back odds
      odds_plac_vec <- rr_to_odds(rr_plac_vec)
      odds_back_vec <- odds_plac_vec * back_add_or
      odds_mono_vec <- odds_plac_vec * mono_add_or
      odds_comb_vec <- odds_plac_vec * comb_add_or
      # transfer odds to rr
      rr_back_vec <- odds_to_rr(odds_back_vec)
      rr_mono_vec <- odds_to_rr(odds_mono_vec)
      rr_comb_vec <- odds_to_rr(odds_comb_vec)
    }

    rr_transform_vec <- rr_transform[sample(1:length(rr_transform), cohorts_max, prob = prob_rr_transform, replace = TRUE)]

  } else {
    rr_comb_vec <- rep(rr_comb, cohorts_max)
    rr_back_vec <- rep(rr_back, cohorts_max)
    rr_mono_vec <- rep(rr_mono, cohorts_max)
    rr_plac_vec <- rep(rr_plac, cohorts_max)
    rr_transform_vec <- rr_transform[sample(1:length(rr_transform), cohorts_max, prob = 1, replace = TRUE)]
  }


  # dummy to indicate trial stop
  trial_stop <- 0

  # dummy for timestamp for first success
  first_success <- -1

  # Variable measuring time since last cohort was added
  last_cohort_time <- 0

  # Initialize res_list
  res_list <- create_cohort_initial(cohorts_start, design_type, arms_per_cohort,
                                    rr_comb_vec, rr_mono_vec, rr_back_vec, rr_plac_vec)

  Total_N_Vector <- NULL

  # initialize platform time
  plat_time <- 0

  # Initialize new Patient vector and vector of exact arrival times
  new_pats <- NULL
  pats_arrival_times <- NULL

  # Initialize empty data frame
  cols <- c("PatID", "ArrivalTime", "Cohort", "Arm", "RespBio", "RespHist", "HistMissing")
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
      new_pats <- c(new_pats, rpois(1, accrual_param))
    }

    if (accrual_type == "exponential") {
      new_pats <- c(new_pats, round(accrual_param ^ (plat_time)))
    }

    new_arrival_times <- sort(runif(new_pats[plat_time])) + plat_time - 1
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
        f <- match.fun(rr_transform_vec[[new_pat_cohort]])
        new_probs <- f(res_list[[new_pat_cohort]]$Arms[[new_pat_arm]]$rr[plat_time])
        draw <- t(stats::rmultinom(1, 1, new_probs))
        resp_bio <- 0
        resp_hist <- 0
        if (draw[1,2] == 1) {
          resp_bio <- resp_bio + 1
        }
        if (draw[1,3] == 1) {
          resp_hist <- resp_hist + 1
        }
        if (draw[1,4] == 1) {
          resp_hist <- resp_hist + 1
          resp_bio <- resp_bio + 1
        }


        # create data of new patient
        single_pat <-
          data.frame(
            PatID       = nrow(df) + 1, # Patient ID is one more than what was there previously
            ArrivalTime = i, # current exact time
            Cohort      = new_pat_cohort,
            Arm         = new_pat_arm,
            RespBio     = resp_bio,
            RespHist    = resp_hist,
            HistMissing = sample(0:1, 1, prob = c(1 - missing_prob, missing_prob)) # sample whether hist endpoint will be missing
          )

        # Add patient to existing data frame
        df <- rbind(df, single_pat)
      }

      # Update the Meta information for every arm/cohort (i.e. how many outcomes are observed)
      res_list <- observe_outcomes(res_list, i, bio_lag, hist_lag)

      # Check whether any analyses should be conducted now
      # i.e. for interim check, if at least n_int patients per cohort have observed the biomarker outcome
      ind_interim <- which(sapply(res_list, function(x) check_int_milestone(x, n_int = n_int)))
      ind_final   <- which(sapply(res_list, function(x) check_fin_milestone(x, n_fin = n_fin)))

      if (length(ind_interim) > 0) {
        ############# conduct interim analysis
        for (j in ind_interim) {

          res_list <-
            make_decision_trial(
              res_list         = res_list,
              which_cohort     = j,
              interim          = TRUE,
              design_type      = design_type,
              arms_per_cohort  = arms_per_cohort,
              analysis_time    = i,
              bio_lag          = bio_lag,
              hist_lag         = hist_lag,
              dataset          = df,
              sharing_type     = sharing_type,
              ...
            )

          # What happens at successful interim
          if (res_list[[j]]$Meta$decision[1] == "GO_SUP") {
            res_list[[j]]$Meta$decision[2] <- "GO_SUP"
            if (first_success == -1) {
              first_success <- plat_time
            }
            rand_list <- create_rand_list(res_list)
          }

          # What happens at unsuccessful interim
          if (res_list[[j]]$Meta$decision[1] == "STOP_FUT") {
            res_list[[j]]$Meta$decision[2] <- "STOP_FUT"
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
              interim          = FALSE,
              design_type      = design_type,
              arms_per_cohort  = arms_per_cohort,
              analysis_time    = i,
              bio_lag          = bio_lag,
              hist_lag         = hist_lag,
              dataset          = df,
              sharing_type     = sharing_type,
              ...
            )

          if (res_list[[j]]$Meta$decision[2] == "GO_SUP") {
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
        res_list[[new_patcohort]]$Meta$decision[2] <- "STOP_SAFETY"
        rand_list <- create_rand_list(res_list)
      }

      # If all cohorts are stopped, stop trial
      if (!any(sapply(res_list, function(x) x$Meta$decision[2]) %in% c("none", "CONTINUE"))) {
        trial_stop <- 1
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

    # make sure that number of arms to add is bounded above by cohorts_max and cohorts_sim
    coh_add <- min(coh_add, cohorts_max - length(res_list), cohorts_sim - sum(sapply(res_list, function(x) coh_left_check(x))))

    # Add new cohorts and create new randomization list immediately
    if (coh_add > 0) {
      for (c in 1:coh_add) {
        res_list <-
          create_cohort_new(
            res_list,
            plat_time,
            design_type,
            arms_per_cohort,
            rr_comb_vec,
            rr_mono_vec,
            rr_back_vec,
            rr_plac_vec
          )
      }

      last_cohort_time <- plat_time
      rand_list <- create_rand_list(res_list)
    }

    # if trial continues, update response rates
    res_list <- update_rr(res_list)

  }

  # Define truth via:
  # a) Risk Difference & Combination Design
  # a1) Combo > SoC (2 arms)
  # a2) Combo > Mono > SoC (3+4 arms)
  # b) Risk Difference & Doses Design
  # b1) Any dose > SoC
  # Use initial RR vectors that were assigned

  truth <- rep(NA, length(res_list))


    for (i in 1:length(res_list)) {
      if (arms_per_cohort == 2) {
        truth[i] <-
          (res_list[[i]]$Arms$Comb$rr[1] > res_list[[i]]$Arms$Plac$rr[1])
      } else if (arms_per_cohort == 3) {

        if (design_type == "combination") {

          (res_list[[i]]$Arms$Comb$rr[1] > res_list[[i]]$Arms$Mono$rr[1]) &
          (res_list[[i]]$Arms$Mono$rr[1] > res_list[[i]]$Arms$Plac$rr[1])

        }

        if (design_type == "doses") {

          (res_list[[i]]$Arms$Comb$rr[1] > res_list[[i]]$Arms$Plac$rr[1]) &
          (res_list[[i]]$Arms$Mono$rr[1] > res_list[[i]]$Arms$Plac$rr[1])

        }

      } else {

        if (design_type == "combination") {

          truth[i] <-
            (res_list[[i]]$Arms$Comb$rr[1] > res_list[[i]]$Arms$Mono$rr[1]) &
            (res_list[[i]]$Arms$Comb$rr[1] > res_list[[i]]$Arms$Back$rr[1]) &
            (res_list[[i]]$Arms$Mono$rr[1] > res_list[[i]]$Arms$Plac$rr[1]) &
            (res_list[[i]]$Arms$Back$rr[1] > res_list[[i]]$Arms$Plac$rr[1])

        }

        if (design_type == "doses") {

          (res_list[[i]]$Arms$Comb$rr[1] > res_list[[i]]$Arms$Plac$rr[1]) &
          (res_list[[i]]$Arms$Mono$rr[1] > res_list[[i]]$Arms$Plac$rr[1]) &
          (res_list[[i]]$Arms$Back$rr[1] > res_list[[i]]$Arms$Plac$rr[1])

        }

      }

  }

  # Get final experimental response rates over time
  rr_comb_final <- sapply(res_list, function(x) x$Arms$Comb$rr)
  rr_plac_final <- sapply(res_list, function(x) x$Arms$Plac$rr)
  rr_mono_final <- sapply(res_list, function(x) x$Arms$Mono$rr)
  rr_back_final <- sapply(res_list, function(x) x$Arms$Back$rr)

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
  cp <- sum(substring(sapply(res_list, function(x) x$Meta$decision[2]), 1, 2) == "GO" &  truth)
  fp <- sum(substring(sapply(res_list, function(x) x$Meta$decision[2]), 1, 2) == "GO" & !truth)
  cn <- sum(substring(sapply(res_list, function(x) x$Meta$decision[2]), 1, 2) == "ST" & !truth)
  fn <- sum(substring(sapply(res_list, function(x) x$Meta$decision[2]), 1, 2) == "ST" &  truth)

  # Prepare return list
  ret <- list(
    Decision               = sapply(res_list, function(x) x$Meta$decision),
    Start_N                = sapply(res_list, function(x) x$Meta$start_n),
    Start_Time             = sapply(res_list, function(x) x$Meta$start_time),
    RR_Comb                = rr_comb_final,
    RR_Mono                = rr_mono_final,
    RR_Back                = rr_back_final,
    RR_Plac                = rr_plac_final,
    N_Cohorts              = length(res_list),
    Final_N_Cohort         = sapply(res_list, function(x) x$Meta$pat_enrolled),
    Final_N_Cohort_Trial   = mean(sapply(res_list, function(x) x$Meta$pat_enrolled)),
    Total_N                = sum(sapply(res_list, function(x) x$Meta$pat_enrolled)),
    Overrunning            = sapply(res_list, function(x) ifelse(x$Meta$decision[1] == "CONTINUE", 0, x$Meta$pat_enrolled - n_int)),
    # does not take into account random stopping
    Overrunning_Trial      = mean(sapply(res_list, function(x) ifelse(x$Meta$decision[1] == "CONTINUE", 0, x$Meta$pat_enrolled - n_int))),
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
    Int_GO                 = sum(sapply(res_list, function(x) x$Meta$sup_interim), na.rm = TRUE),
    Int_STOP               = sum(sapply(res_list, function(x) x$Meta$fut_interim), na.rm = TRUE),
    Safety_STOP            = sum(sapply(res_list, function(x) (x$Meta$decision[2] == "STOP_SAFETY")), na.rm = TRUE),
    Int_GO_Trial           = sum(sapply(res_list, function(x) x$Meta$sup_interim), na.rm = TRUE) / length(res_list),
    Int_STOP_Trial         = sum(sapply(res_list, function(x) x$Meta$fut_interim), na.rm = TRUE) / length(res_list),
    Safety_STOP_Trial      = sum(sapply(res_list, function(x) (x$Meta$decision[2] == "STOP_SAFETY")), na.rm = TRUE) / length(res_list)
  )

  if (stage_data) {
    ret <- list(Trial_Overview = ret, Stage_Data = res_list, Pat_Data = df)
  } else {
    ret <- list(Trial_Overview = ret)
  }

  # Whatif Analysis

  whatif_dec <- rep(NA, length(res_list))

  for (i in 1:length(whatif_dec)) {

    res_list <-
      make_decision_trial(
        res_list         = res_list,
        which_cohort     = i,
        interim          = FALSE,
        design_type      = design_type,
        arms_per_cohort  = arms_per_cohort,
        analysis_time    = Inf,
        bio_lag          = bio_lag,
        hist_lag         = hist_lag,
        dataset          = df,
        hist_missing     = FALSE,
        sharing_type     = sharing_type,
        ...
      )

    whatif_dec[i] <- res_list[[i]]$Meta$decision[2]

  }

  ret$Trial_Overview <-
    c(
      Whatif = list(whatif_dec),
      Whatif_Agreement_Trial = mean(whatif_dec == ret$Trial_Overview$Decision[2,]),
      ret$Trial_Overview
  )

  return(ret)

}
