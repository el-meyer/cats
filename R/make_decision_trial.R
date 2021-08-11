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
#' @param Bayes_Sup        List of matrices with rows corresponding to number of multiple Bayesian posterior two-arm combination criteria for superiority
#'
#' @param Bayes_Fut        List of matrices with rows corresponding to number of multiple Bayesian posterior two-arm combination criteria for futility
#'
#' @param P_Sup            List with sublists corresponding to number of multiple frequentist test-based combination criteria for superiority
#'
#' @param P_Fut            List with sublists corresponding to number of multiple frequentist test-based combination criteria for futility
#'
#' @param interim          Is the analysis conducted an interim or a final analysis?
#'
#' @param analysis_time    Platform Time of Analysis
#'
#' @param dataset          Dataset to be used for analysis
#'
#' @param hist_miss        Whether or not to exclude missing histology data
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
#' cols <- c("PatID", "ArrivalTime", "Cohort", "Arm", "RespBio", "RespHist", "HistMissing")
#' df <- matrix(nrow = 100, ncol = length(cols))
#' colnames(df) <- cols
#' df <- as.data.frame(df)
#' df$PatID <- 1:100
#' df$ArrivalTime <- sort(runif(100, min = 0, max = 5))
#' df$Cohort <- sample(1:2, 100, replace = TRUE)
#' df$Arm <- sample(c("Combo", "Plac"), 100, replace = TRUE)
#' df$RespBio <- sample(0:1, 100, replace = TRUE)
#' df$RespHist <- sample(0:1, 100, replace = TRUE)
#' df$HistMissing <- sample(0:1, 100, replace = TRUE, prob = c(0.95, 0.05))
#'
#' # Comparison Combo vs Mono
#' Bayes_Sup1 <- matrix(nrow = 3, ncol = 2)
#' Bayes_Sup1[1,] <- c(0.00, 0.95)
#' Bayes_Sup1[2,] <- c(0.10, 0.80)
#' Bayes_Sup1[3,] <- c(0.15, 0.50)
#' # Comparison Combo vs Backbone
#' Bayes_Sup2 <- matrix(nrow = 3, ncol = 2)
#' Bayes_Sup2[1,] <- c(0.00, 0.95)
#' Bayes_Sup2[2,] <- c(NA, NA)
#' Bayes_Sup2[3,] <- c(NA, NA)
#' # Comparison Mono vs Placebo
#' Bayes_Sup3 <- matrix(nrow = 3, ncol = 2)
#' Bayes_Sup3[1,] <- c(0.00, 0.95)
#' Bayes_Sup3[2,] <- c(0.10, 0.80)
#' Bayes_Sup3[3,] <- c(NA, NA)
#' #' # Comparison Backbone vs Placebo
#' Bayes_Sup4 <- matrix(nrow = 3, ncol = 2)
#' Bayes_Sup4[1,] <- c(0.00, 0.95)
#' Bayes_Sup4[2,] <- c(0.10, 0.80)
#' Bayes_Sup4[3,] <- c(NA, NA)
#' Bayes_Sup <- list(list(Bayes_Sup1, Bayes_Sup2, Bayes_Sup3, Bayes_Sup4),
#'                   list(Bayes_Sup1, Bayes_Sup2, Bayes_Sup3, Bayes_Sup4))
#'
#' sharing_type <- "all"
#' interim <- TRUE
#' which_cohort <- 1
#'
#' # DO NOT RUN
#' make_decision_trial(
#' res_list = res_list, which_cohort = which_cohort,
#' interim = interim, missing_prob = missing_prob,
#' Bayes_Sup = Bayes_Sup, sharing_type = sharing_type,
#' seed_missing = seed_missing,
#' )
#'
#' # Multiple decision rules
#'
#' # Vergleich Combo vs Mono
#' Bayes_Fut1 <- matrix(nrow = 1, ncol = 2)
#' Bayes_Fut1[1,] <- c(NA, NA)
#' # Vergleich Combo vs Backbone
#' Bayes_Fut2 <- matrix(nrow = 1, ncol = 2)
#' Bayes_Fut2[1,] <- c(NA, NA)
#' # Vergleich Mono vs Placebo
#' Bayes_Fut3 <- matrix(nrow = 1, ncol = 2)
#' Bayes_Fut3[1,] <- c(0.00, 0.60)
#' Bayes_Fut4 <- matrix(nrow = 1, ncol = 2)
#' Bayes_Fut4[1,] <- c(0.00, 0.60)
#' Bayes_Fut <- list(list(Bayes_Fut1, Bayes_Fut2, Bayes_Fut3, Bayes_Fut4),
#'                   list(Bayes_Fut1, Bayes_Fut2, Bayes_Fut3, Bayes_Fut4))
#'
#' # Comparison Combo vs Mono
#' P_Sup1 <- list(list(
#' testfun = function(x) stats::prop.test(x, alternative = "less", correct = FALSE),
#' p_sup = 0.025, p_adj = "B"))
#' # Comparison Combo vs Backbone
#' P_Sup2 <- list(list(
#' testfun = function(x) stats::prop.test(x, alternative = "less", correct = FALSE),
#' p_sup = 0.025, p_adj = "B"))
#' # Comparison Mono vs Placebo
#' P_Sup3 <- list(list(
#' testfun = function(x) stats::prop.test(x, alternative = "less", correct = FALSE),
#' p_sup = 0.050, p_adj = "B"))
#' P_Sup4 <- list(list(
#' testfun = function(x) stats::prop.test(x, alternative = "less", correct = FALSE),
#' p_sup = 0.050, p_adj = "B"))
#' P_Sup <- list(list(P_Sup1, P_Sup2, P_Sup3, P_Sup4),
#'               list(P_Sup1, P_Sup2, P_Sup3, P_Sup4))
#'
#' # Comparison Combo vs Mono
#' P_Fut1 <- list(list(
#' testfun = function(x) stats::prop.test(x, alternative = "less", correct = FALSE),
#' p_fut = 0.5, p_adj = "none"))
#' # Comparison Combo vs Backbone
#' P_Fut2 <- list(list(
#' testfun = function(x) stats::prop.test(x, alternative = "less", correct = FALSE),
#' p_fut = 0.5, p_adj = "none"))
#' # Comparison Mono vs Placebo
#' P_Fut3 <- list(list(
#' testfun = function(x) stats::prop.test(x, alternative = "less", correct = FALSE),
#' p_fut = 0.5, p_adj = "none"))
#' # Comparison Backbone Placebo
#' P_Fut4 <- list(list(
#' testfun = function(x) stats::prop.test(x, alternative = "less", correct = FALSE),
#' p_fut = 0.5, p_adj = "none"))
#' P_Fut <- list(list(P_Fut1, P_Fut2, P_Fut3, P_Fut4),
#'               list(P_Fut1, P_Fut2, P_Fut3, P_Fut4))
#'
#' # DO NOT RUN
#' make_decision_trial(res_list = res_list, which_cohort = which_cohort, interim = interim,
#' Bayes_Sup = Bayes_Sup, sharing_type = sharing_type,
#' Bayes_Fut = Bayes_Fut, Bayes_SA_Sup = Bayes_SA_Sup, Bayes_SA_Fut = Bayes_SA_Fut, P_Sup = P_Sup,
#' P_Fut = P_Fut, Est_Sup_Fut = Est_Sup_Fut, CI_Sup_Fut = CI_Sup_Fut
#' )
#'
#' @export
make_decision_trial <- function(res_list, which_cohort, Bayes_Sup = NULL, Bayes_Fut = NULL,
                                w = 0.5, P_Sup = NULL, P_Fut = NULL, interim, beta_prior = 0.5,
                                analysis_time, dataset, hist_miss = TRUE, ...) {

  # todo: allow Estimates and CIs to be single arm rules for combo, mono and back
  # (atm only implemented as two-arm decision rules)

  ########## Parameter Explanation and Variable and Function Initiation ##########

  # Bayes_Sup = Matrix with rows corresponding to number of multiple Bayesian two-arm posterior combination criteria for superiority
  # First element is always the delta, second element is the GO-threshold
  # The rule is: If P(P_e > P_c + delta) > GO-threshold, then GO
  # The second rule is: If promising_threshold < P(P_e > P_c + delta) < GO-threshold, then declare promising

  # Bayes_Fut = Matrix with rows corresponding to number of multiple Bayesian two-arm posterior combination criteria for futility
  # First element is always the delta and the second element is the STOP-threshold
  # The rule is: If P(P_e > P_c + delta) < STOP-threshold, then STOP

  # P_Sup = List with sublists corresponding to number of multiple frequentist test-based combination criteria for superiority
  # Every sublist has the same structure:
  # The first element is the Test Function. This needs to be a function out of a list of allowed functions, or any
  # custom function that takes a 2x2 contingency table as input
  # The second element is the Superiority Threshold (i.e. GO, if p < p_sup)
  # The third element contains the p-value multiplicity adjustment method at every decision

  # P_Fut = List with sublists corresponding to number of multiple frequentist test-based combination criteria for futility
  # Every sublist has the same structure:
  # The first element is the Test Function. This needs to be a function out of a list of allowed functions, or any
  # custom function that takes a 2x2 contingency table as input
  # The second element is the Superiority Threshold (i.e. STOP, if p > p_fut)

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

  "%>%" <- dplyr::"%>%"

  # Create list for whether futility/superiority decisions were made for every single rule for all three comparisons
  # All list elements will be matrices, with rows corresponding to multiple decision rules and columns corresponding to the
  # three comparisons (j=1: Combo vs Mono, j=2: Combo vs Back, j=3: Mono vs Placebo, j=4: Back vs Placebo)
  superiority_reached <- list()
  futility_reached <- list()

  sup_reached <- fut_reached <- NULL

  # Get responders and sample sizes
  resp_bio_tot <- resp_hist_tot <- n_tot <- NULL

  # Combo
  if (interim) {

    resp_bio <-
      dataset %>%
      dplyr::filter(
        Cohort == which_cohort,
        Arm == "Comb",
        ArrivalTime < analysis_time - bio_lag
      ) %>%
      dplyr::pull(RespBio)

    resp_bio_tot[1] <- sum(resp_bio)
    n_tot[1] <- length(resp_bio)

  } else {

    resp_hist <-
      dataset %>%
      dplyr::filter(
        Cohort == which_cohort,
        Arm == "Comb",
        ArrivalTime < analysis_time - hist_lag # should be all
      ) %>%
        {if (hist_miss) dplyr::filter(., HistMissing != 1) else .} %>%
      dplyr::pull(RespHist)

    resp_hist_tot[1] <- sum(resp_hist)
    n_tot[1] <- length(resp_hist)

  }

  # Plac
  # Indicator of position is "arms_per_cohort"
  if (interim) {

    # if sharing all data, time frame to take under consideration is simply all data from all cohorts
    if (sharing_type == "all") {
      resp_bio <-
        dataset %>%
        dplyr::filter(
          Arm == "Plac",
          ArrivalTime < analysis_time - bio_lag
        ) %>%
        dplyr::pull(RespBio)

      resp_bio_tot[arms_per_cohort] <- sum(resp_bio)
      n_tot[arms_per_cohort] <- length(resp_bio)
    }

    # if sharing no data, take only data from current Arm
    if (sharing_type == "cohort") {
      resp_bio <-
        dataset %>%
        dplyr::filter(
          Cohort == which_cohort,
          Arm == "Plac",
          ArrivalTime < analysis_time - bio_lag
        ) %>%
        dplyr::pull(RespBio)

      resp_bio_tot[arms_per_cohort] <- sum(resp_bio)
      n_tot[arms_per_cohort] <- length(resp_bio)
    }

    # if sharing concurrent data, take only data from certain time period
    if (sharing_type == "concurrent") {
      resp_bio <-
        dataset %>%
        dplyr::filter(
          Arm == "Plac",
          ArrivalTime < analysis_time - bio_lag,
          ArrivalTime > res_list[[which_cohort]]$Meta$start_time
        ) %>%
        dplyr::pull(RespBio)

      resp_bio_tot[arms_per_cohort] <- sum(resp_bio)
      n_tot[arms_per_cohort] <- length(resp_bio)
    }

    # if dynamic borrowing, create two datasets and combine
    if (sharing_type == "dynamic") {
      resp_bio_cohort <-
        dataset %>%
        dplyr::filter(
          Cohort == which_cohort,
          Arm == "Plac",
          ArrivalTime < analysis_time - bio_lag
        ) %>%
        dplyr::pull(RespBio)

      resp_bio_external <-
        dataset %>%
        dplyr::filter(
          Cohort != which_cohort,
          Arm == "Plac",
          ArrivalTime < analysis_time - bio_lag
        ) %>%
        dplyr::pull(RespBio)

      # compute responders from cohort
      suc_bio_c <- sum(resp_bio_cohort)
      N_c <- length(resp_bio_cohort)

      # compute historical responders
      suc_bio_h <- sum(resp_bio_external)
      N_h <- length(resp_bio_external)

      w1_bio <- w * beta(suc_bio_c + suc_bio_h + beta_prior, N_h + N_c - suc_bio_c - suc_bio_h + beta_prior) /
        beta(suc_bio_h + beta_prior, N_h - suc_bio_h + beta_prior)
      w2_bio <- (1 - w) * beta(suc_bio_c + beta_prior, N_c - suc_bio_c + beta_prior) / beta(beta_prior, beta_prior)
      w1n_bio <- w1_bio / (w1_bio + w2_bio)
      w2n_bio <- w2_bio / (w1_bio + w2_bio)

      m_a_bio <- w1n_bio * (suc_bio_c + suc_bio_h + beta_prior) + w2n_bio * (suc_bio_c + beta_prior)
      m_b_bio <- w1n_bio * (N_h + N_c - suc_bio_c - suc_bio_h + beta_prior) + w2n_bio * (N_c - suc_bio_c + beta_prior)

      n_tot[arms_per_cohort] <- round(m_a_bio + m_b_bio)
      resp_bio_tot[arms_per_cohort] <- round(m_a_bio)
    }

  } else {

    # if sharing all data, time frame to take under consideration is simply all data from all cohorts
    if (sharing_type == "all") {
      resp_hist <-
        dataset %>%
        dplyr::filter(
          Arm == "Plac",
          ArrivalTime < analysis_time - hist_lag # should be all
        ) %>%
        {if (hist_miss) dplyr::filter(., HistMissing != 1) else .} %>%
        dplyr::pull(RespHist)

      resp_hist_tot[arms_per_cohort] <- sum(resp_hist)
      n_tot[arms_per_cohort] <- length(resp_hist)
    }

    # if sharing no data, take only data from current Arm
    if (sharing_type == "cohort") {
      resp_hist <-
        dataset %>%
        dplyr::filter(
          Cohort == which_cohort,
          Arm == "Plac",
          ArrivalTime < analysis_time - hist_lag # should be all
        ) %>%
        {if (hist_miss) dplyr::filter(., HistMissing != 1) else .} %>%
        dplyr::pull(RespHist)

      resp_hist_tot[arms_per_cohort] <- sum(resp_hist)
      n_tot[arms_per_cohort] <- length(resp_hist)
    }

    # if sharing concurrent data, take only data from certain time period
    if (sharing_type == "concurrent") {
      resp_hist <-
        dataset %>%
        dplyr::filter(
          Arm == "Plac",
          ArrivalTime < analysis_time - hist_lag,
          ArrivalTime > res_list[[which_cohort]]$Meta$start_time,
        ) %>%
        {if (hist_miss) dplyr::filter(., HistMissing != 1) else .} %>%
        dplyr::pull(RespHist)

      resp_hist_tot[arms_per_cohort] <- sum(resp_hist)
      n_tot[arms_per_cohort] <- length(resp_hist)
    }

    # if dynamic borrowing, create two datasets and combine
    if (sharing_type == "dynamic") {
      resp_hist_cohort <-
        dataset %>%
        dplyr::filter(
          Cohort == which_cohort,
          Arm == "Plac",
          ArrivalTime < analysis_time - hist_lag # should be all
        ) %>%
        {if (hist_miss) dplyr::filter(., HistMissing != 1) else .} %>%
        dplyr::pull(RespHist)

      resp_hist_external <-
        dataset %>%
        dplyr::filter(
          Cohort != which_cohort,
          Arm == "Plac",
          ArrivalTime < analysis_time - hist_lag # should be all
        ) %>%
        {if (hist_miss) dplyr::filter(., HistMissing != 1) else .} %>%
        dplyr::pull(RespHist)

      # compute responders from cohort
      suc_bio_c <- sum(resp_bio_cohort)
      N_c <- length(resp_bio_cohort)

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

      n_tot[arms_per_cohort] <- round(m_a_hist + m_b_hist)
      resp_hist_tot[arms_per_cohort] <- round(m_a_hist)
    }

  }

  # if arms per cohort is 3, add mono/medium

  if (arms_per_cohort == 3) {

    # Mono
    if (interim) {

    resp_bio <-
        dataset %>%
        dplyr::filter(
          Cohort == which_cohort,
          Arm == "Mono",
          ArrivalTime < analysis_time - bio_lag
        ) %>%
        dplyr::pull(RespBio)

      resp_bio_tot[2] <- sum(resp_bio)
      n_tot[2] <- length(resp_bio)

    } else {

      resp_hist <-
        dataset %>%
        dplyr::filter(
          Cohort == which_cohort,
          Arm == "Mono",
          ArrivalTime < analysis_time - hist_lag # should be all
        ) %>%
        {if (hist_miss) dplyr::filter(., HistMissing != 1) else .} %>%
        dplyr::pull(RespHist)

      resp_hist_tot[2] <- sum(resp_hist)
      n_tot[2] <- length(resp_hist)

    }

  }

  # if arms per cohort is 4, add mono/medium and backbone/low

  if (arms_per_cohort == 4) {

    # Mono
    if (interim) {

      resp_bio <-
        dataset %>%
        dplyr::filter(
          Cohort == which_cohort,
          Arm == "Mono",
          ArrivalTime < analysis_time - bio_lag
        ) %>%
        dplyr::pull(RespBio)

      resp_bio_tot[2] <- sum(resp_bio)
      n_tot[2] <- length(resp_bio)

    } else {

      resp_hist <-
        dataset %>%
        dplyr::filter(
          Cohort == which_cohort,
          Arm == "Mono",
          ArrivalTime < analysis_time - hist_lag # should be all
        ) %>%
        {if (hist_miss) dplyr::filter(., HistMissing != 1) else .} %>%
        dplyr::pull(RespHist)

      resp_hist_tot[2] <- sum(resp_hist)
      n_tot[2] <- length(resp_hist)

    }

    # Backbone
    # differentiate between trial design:
    # In dose design, no sharing
    # In combination design, do sharing

    if (design_type == "combination") {

      if (interim) {

        # if sharing all data, time frame to take under consideration is simply all data from all cohorts
        if (sharing_type == "all") {
          resp_bio <-
            dataset %>%
            dplyr::filter(
              Arm == "Back",
              ArrivalTime < analysis_time - bio_lag
            ) %>%
            dplyr::pull(RespBio)

          resp_bio_tot[arms_per_cohort] <- sum(resp_bio)
          n_tot[arms_per_cohort] <- length(resp_bio)
        }

        # if sharing no data, take only data from current Arm
        if (sharing_type == "cohort") {
          resp_bio <-
            dataset %>%
            dplyr::filter(
              Cohort == which_cohort,
              Arm == "Back",
              ArrivalTime < analysis_time - bio_lag
            ) %>%
            dplyr::pull(RespBio)

          resp_bio_tot[arms_per_cohort] <- sum(resp_bio)
          n_tot[arms_per_cohort] <- length(resp_bio)
        }

        # if sharing concurrent data, take only data from certain time period
        if (sharing_type == "concurrent") {
          resp_bio <-
            dataset %>%
            dplyr::filter(
              Arm == "Back",
              ArrivalTime < analysis_time - bio_lag,
              ArrivalTime > res_list[[which_cohort]]$Meta$start_time
            ) %>%
            dplyr::pull(RespBio)

          resp_bio_tot[arms_per_cohort] <- sum(resp_bio)
          n_tot[arms_per_cohort] <- length(resp_bio)
        }

        # if dynamic borrowing, create two datasets and combine
        if (sharing_type == "dynamic") {
          resp_bio_cohort <-
            dataset %>%
            dplyr::filter(
              Cohort == which_cohort,
              Arm == "Back",
              ArrivalTime < analysis_time - bio_lag
            ) %>%
            dplyr::pull(RespBio)

          resp_bio_external <-
            dataset %>%
            dplyr::filter(
              Cohort != which_cohort,
              Arm == "Back",
              ArrivalTime < analysis_time - bio_lag
            ) %>%
            dplyr::pull(RespBio)

          # compute responders from cohort
          suc_bio_c <- sum(resp_bio_cohort)
          N_c <- length(resp_bio_cohort)

          # compute historical responders
          suc_bio_h <- sum(resp_bio_external)
          N_h <- length(resp_bio_external)

          w1_bio <- w * beta(suc_bio_c + suc_bio_h + beta_prior, N_h + N_c - suc_bio_c - suc_bio_h + beta_prior) /
            beta(suc_bio_h + beta_prior, N_h - suc_bio_h + beta_prior)
          w2_bio <- (1 - w) * beta(suc_bio_c + beta_prior, N_c - suc_bio_c + beta_prior) / beta(beta_prior, beta_prior)
          w1n_bio <- w1_bio / (w1_bio + w2_bio)
          w2n_bio <- w2_bio / (w1_bio + w2_bio)

          m_a_bio <- w1n_bio * (suc_bio_c + suc_bio_h + beta_prior) + w2n_bio * (suc_bio_c + beta_prior)
          m_b_bio <- w1n_bio * (N_h + N_c - suc_bio_c - suc_bio_h + beta_prior) + w2n_bio * (N_c - suc_bio_c + beta_prior)

          n_tot[arms_per_cohort] <- round(m_a_bio + m_b_bio)
          resp_bio_tot[arms_per_cohort] <- round(m_a_bio)
        }

      } else {

        # if sharing all data, time frame to take under consideration is simply all data from all cohorts
        if (sharing_type == "all") {
          resp_hist <-
            dataset %>%
            dplyr::filter(
              Arm == "Back",
              ArrivalTime < analysis_time - hist_lag # should be all
            ) %>%
            {if (hist_miss) dplyr::filter(., HistMissing != 1) else .} %>%
            dplyr::pull(RespHist)

          resp_hist_tot[arms_per_cohort] <- sum(resp_hist)
          n_tot[arms_per_cohort] <- length(resp_hist)
        }

        # if sharing no data, take only data from current Arm
        if (sharing_type == "cohort") {
          resp_hist <-
            dataset %>%
            dplyr::filter(
              Cohort == which_cohort,
              Arm == "Back",
              ArrivalTime < analysis_time - hist_lag # should be all
            ) %>%
            {if (hist_miss) dplyr::filter(., HistMissing != 1) else .} %>%
            dplyr::pull(RespHist)

          resp_hist_tot[arms_per_cohort] <- sum(resp_hist)
          n_tot[arms_per_cohort] <- length(resp_hist)
        }

        # if sharing concurrent data, take only data from certain time period
        if (sharing_type == "concurrent") {
          resp_hist <-
            dataset %>%
            dplyr::filter(
              Arm == "Back",
              ArrivalTime < analysis_time - hist_lag,
              ArrivalTime > res_list[[which_cohort]]$Meta$start_time,
            ) %>%
            {if (hist_miss) dplyr::filter(., HistMissing != 1) else .} %>%
            dplyr::pull(RespHist)

          resp_hist_tot[arms_per_cohort] <- sum(resp_hist)
          n_tot[arms_per_cohort] <- length(resp_hist)
        }

        # if dynamic borrowing, create two datasets and combine
        if (sharing_type == "dynamic") {
          resp_hist_cohort <-
            dataset %>%
            dplyr::filter(
              Cohort == which_cohort,
              Arm == "Back",
              ArrivalTime < analysis_time - hist_lag # should be all
            ) %>%
            {if (hist_miss) dplyr::filter(., HistMissing != 1) else .} %>%
            dplyr::pull(RespHist)

          resp_hist_external <-
            dataset %>%
            dplyr::filter(
              Cohort != which_cohort,
              Arm == "Back",
              ArrivalTime < analysis_time - hist_lag,
              HistMissing != 1
            ) %>%
            dplyr::pull(RespHist)

          # compute responders from cohort
          suc_bio_c <- sum(resp_bio_cohort)
          N_c <- length(resp_bio_cohort)

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

          n_tot[arms_per_cohort] <- round(m_a_hist + m_b_hist)
          resp_hist_tot[arms_per_cohort] <- round(m_a_hist)
        }

      }

    }

    if (design_type == "doses") {

      if (interim) {

        resp_bio <-
          dataset %>%
          dplyr::filter(
            Cohort == which_cohort,
            Arm == "Back",
            ArrivalTime < analysis_time - bio_lag
          ) %>%
          dplyr::pull(RespBio)

        resp_bio_tot[3] <- sum(resp_bio)
        n_tot[3] <- length(resp_bio)

      } else {

        resp_hist <-
          dataset %>%
          dplyr::filter(
            Cohort == which_cohort,
            Arm == "Back",
            ArrivalTime < analysis_time - hist_lag # should be all
          ) %>%
          {if (hist_miss) dplyr::filter(., HistMissing != 1) else .} %>%
          dplyr::pull(RespHist)

        resp_hist_tot[3] <- sum(resp_hist)
        n_tot[3] <- length(resp_hist)

      }

    }

  }


  if (interim) {
    resp_tot <- resp_bio_tot
  } else {
    resp_tot <- resp_hist_tot
  }

  # define number of comparisons and "contrast"
  num_comp <- NULL
  if (design_type == "doses") {
    num_comp <- arms_per_cohort - 1
  }
  if (design_type == "combination") {
    if (arms_per_cohort == 2) {
      num_comp <- 1
    } else if (arms_per_cohort == 3) {
      num_comp <- 2
    } else {
      num_comp <- 4
    }
  }
  comp <- matrix(ncol = 2, nrow = num_comp)
  colnames(comp) <- c("c1", "c2")
  if (nrow(comp) == 1) {
    comp[,1] <- 1
    comp[,2] <- 2
  } else if (nrow(comp) == 2) {
    if (design_type == "doses") {
      comp[,1] <- c(1,2)
      comp[,2] <- c(3,3)
    } else {
      comp[,1] <- c(1,2)
      comp[,2] <- c(2,3)
    }
  } else if (nrow(comp) == 3) {
    comp[,1] <- c(1,2,3)
    comp[,2] <- c(4,4,4)
  } else {
    comp[,1] <- c(1,1,2,3)
    comp[,2] <- c(2,3,4,4)
  }


  ########## Bayesian Two-Arm Superiority Criteria ###############

  if (!is.null(Bayes_Sup)) {
    # Get the decision rule for interim/final
    if (interim) {
      Bayes_Sup <- Bayes_Sup[[1]]
    } else {
      Bayes_Sup <- Bayes_Sup[[2]]
    }

    # Firstly initiate as many vectors as will maximum be needed (which is determined by the most complex decision rule)
    for (k in 1:max(sapply(Bayes_Sup, function(x) nrow(x)))) {
      assign(paste0("prob_sup_vec", k), rep(NA, num_comp))
    }

    assign(paste0("sup_reached"), matrix(NA, nrow = max(sapply(Bayes_Sup, function(x) nrow(x))), ncol = num_comp))

    # For each of the comparisons separately, evaluate the applying decision rules
    # j=1: Combo vs Mono, j=2: Combo vs Back, j=3: Mono vs Placebo

    for (i in 1:max(sapply(Bayes_Sup, function(x) nrow(x)))) {
      for (j in 1:num_comp) {
        # Get which type of comparison should be conducted
        c1 <- comp[j,1]
        c2 <- comp[j,2]

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


      # Write results in res_list
      eval(parse(text = paste(paste0("res_list[[which_cohort]]$Meta$prob_sup", i, "<-",
                                     "c(res_list[[which_cohort]]$Meta$prob_sup", i, ",prob_sup_vec", i, ")"))))
    }

    # Add decisions from all decision rules for each arm together
    superiority_reached <- c(superiority_reached, list(sup_reached))

  }


  ########## Bayesian Two-Arm Futility Criteria ##########

  if (!is.null(Bayes_Fut)) {
    # Get the decision rule for interim/final
    if (interim) {
      Bayes_Fut <- Bayes_Fut[[1]]
    } else {
      Bayes_Fut <- Bayes_Fut[[2]]
    }

    # Firstly initiate as many vectors as will maximum be needed (which is determined by the most complex decision rule)
    for (k in 1:max(sapply(Bayes_Fut, function(x) nrow(x)))) {
      assign(paste0("prob_fut_vec", k), rep(NA, num_comp))
    }

    assign(paste0("fut_reached"), matrix(NA, nrow = max(sapply(Bayes_Fut, function(x) nrow(x))), ncol = num_comp))

    # For each of the comparisons separately, evaluate the applying decision rules
    # j=1: Combo vs Mono, j=2: Combo vs Back, j=3: Mono vs Placebo

    for (i in 1:max(sapply(Bayes_Fut, function(x) nrow(x)))) {
      for (j in 1:num_comp) {
        # Get which type of comparison should be conducted
        c1 <- comp[j,1]
        c2 <- comp[j,2]

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
      # Write results in res_list
      eval(parse(text = paste(paste0("res_list[[which_cohort]]$Meta$prob_fut", i, "<-",
                                     "c(res_list[[which_cohort]]$Meta$prob_fut", i, ",prob_fut_vec", i, ")"))))
    }

    # Add decisions from all decision rules for each arm together
    futility_reached <- c(futility_reached, list(fut_reached))

  }

  ########## P-Value Superiority Criteria ##########

  if (!is.null(P_Sup)) {
    # Get the decision rule for interim/final
    if (interim) {
      P_Sup <- P_Sup[[1]]
    } else {
      P_Sup <- P_Sup[[2]]
    }

    # Firstly initiate as many vectors as will maximum be needed (which is determined by the most complex decision rule)
    for (k in 1:max(sapply(P_Sup, function(x) length(x)))) {
      assign(paste0("p_values", k), rep(NA, num_comp))
      assign(paste0("p_adjusted", k), rep(NA, num_comp))
    }

    assign(paste0("sup_reached"), matrix(NA, nrow = max(sapply(P_Sup, function(x) length(x))), ncol = num_comp))

    # For each of the comparisons separately, evaluate the applying decision rules

    for (i in 1:max(sapply(P_Sup, function(x) length(x)))) {
      for (j in 1:num_comp) {
        # Get which type of comparison should be conducted
        c1 <- comp[j,1]
        c2 <- comp[j,2]

        # Evaluate decision rule, but if non-existant, give NA
        if (is.na(P_Sup[[j]][[i]]$p_sup)) {
          eval(parse(text = paste(paste0("p_adjusted", i, "[", j, "]"), "<- NA")))
        } else {
          f <- match.fun(P_Sup[[j]][[i]]$testfun)

          crosstab <-  matrix(c(resp_tot[c2], n_tot[c2] - resp_tot[c2], resp_tot[c1], n_tot[c1] - resp_tot[c1]),
                              nrow = 2, ncol = 2, byrow = TRUE)
          eval(parse(text = paste(paste0("p_values", i, "[", j, "]"), "<-", f(crosstab)$p.value)))

          # correct for multiplicity
          if (P_Sup[[j]][[i]]$p_adj == "none") {
            eval(parse(text = paste(paste0("p_adjusted", i, "[", j, "]"), "<-", paste0("p_values", i, "[", j, "]"))))
          }
          if (P_Sup[[j]][[i]]$p_adj == "B") {
            eval(parse(text = paste(paste0("p_adjusted", i, "[", j, "]"), "<-", paste0("min(1, 2*p_values", i, "[", j, "])"))))
          }

        }

        # Check whether superiority reached
        eval(parse(text = paste(paste0("sup_reached[", i, ",", j, "]"), "<-",
                                get(paste0("p_adjusted", i))[j] < P_Sup[[j]][[i]]$p_sup)))

      }
      # Write results in res_list
      eval(parse(text = paste(paste0("res_list[[which_cohort]]$Meta$p_values_sup", i, "<-",
                                     "c(res_list[[which_cohort]]$Meta$p_values_sup", i, ",p_adjusted", i, ")"))))
    }

    # Add decisions from all decision rules for each arm together
    superiority_reached <- c(superiority_reached, list(sup_reached))

  }


  ########## P-Value Futility Criteria ##########

  if (!is.null(P_Fut)) {
    # Get the decision rule for interim/final
    if (interim) {
      P_Fut <- P_Fut[[1]]
    } else {
      P_Fut <- P_Fut[[2]]
    }

    # Firstly initiate as many vectors as will maximum be needed (which is determined by the most complex decision rule)
    for (k in 1:max(sapply(P_Fut, function(x) length(x)))) {
      assign(paste0("p_values", k), rep(NA, num_comp))
      assign(paste0("p_adjusted", k), rep(NA, num_comp))
    }

    assign(paste0("fut_reached"), matrix(NA, nrow = max(sapply(P_Fut, function(x) length(x))), ncol = num_comp))

    # For each of the comparisons separately, evaluate the applying decision rules
    # j=1: Combo vs Mono, j=2: Combo vs Back, j=3: Mono vs Placebo

    for (i in 1:max(sapply(P_Fut, function(x) length(x)))) {
      for (j in 1:num_comp) {
        # Get which type of comparison should be conducted
        c1 <- comp[j,1]
        c2 <- comp[j,2]

        # Evaluate decision rule, but if non-existant, give NA
        if (is.na(P_Fut[[j]][[i]]$p_fut)) {
          eval(parse(text = paste(paste0("p_adjusted", i, "[", j, "]"), "<- NA")))
        } else {
          f <- match.fun(P_Fut[[j]][[i]]$testfun)

          crosstab <-  matrix(c(resp_tot[c2], n_tot[c2] - resp_tot[c2], resp_tot[c1], n_tot[c1] - resp_tot[c1]),
                              nrow = 2, ncol = 2, byrow = TRUE)
          eval(parse(text = paste(paste0("p_values", i, "[", j, "]"), "<-", f(crosstab)$p.value)))

          # correct for multiplicity
          if (P_Fut[[j]][[i]]$p_adj == "none") {
            eval(parse(text = paste(paste0("p_adjusted", i, "[", j, "]"), "<-", paste0("p_values", i, "[", j, "]"))))
          }
          if (P_Fut[[j]][[i]]$p_adj == "B") {
            eval(parse(text = paste(paste0("p_adjusted", i, "[", j, "]"), "<-", paste0("min(1, 2*p_values", i, "[", j, "])"))))
          }

        }

        # Check whether futility reached
        eval(parse(text = paste(paste0("fut_reached[", i, ",", j, "]"), "<-",
                                get(paste0("p_adjusted", i))[j] >= P_Fut[[j]][[i]]$p_fut)))
      }
      # Write results in res_list
      eval(parse(text = paste(paste0("res_list[[which_cohort]]$Meta$p_values_fut", i, "<-",
                                     "c(res_list[[which_cohort]]$Meta$p_values_fut", i, ",p_adjusted", i, ")"))))
    }

    # Add decisions from all decision rules for each arm together
    futility_reached <- c(futility_reached, list(fut_reached))

  }


  ########## Combine Results And Return List ##########

  # Combine results: Choose superiority, futility and evaluating according to rules specified
  # Look at futility first (also random safety stop probability)! Then at superiority, then at promising.

  # Helper Function
  check_fun <- function(x, type, design_type) {
    if (length(x) == 0 | all(sapply(x, function(y) all(is.na(y))))) {
      ret <- FALSE
    } else {
      if (type != "fut" & design_type == "combination") {
        ret <- all(sapply(x, function(y) all(y, na.rm = TRUE)), na.rm = TRUE)
      } else {
        ret <- any(sapply(x, function(y) any(y, na.rm = TRUE)), na.rm = TRUE)
      }
    }
    return(ret)
  }

  # Get final results, depending on whether interim or final
  # Currently all single decisions need to be TRUE in order for the total to be TRUE
  if (interim) {
    res_list[[which_cohort]]$Meta$sup_interim_list <- superiority_reached
    res_list[[which_cohort]]$Meta$fut_interim_list <- futility_reached
    res_list[[which_cohort]]$Meta$sup_interim <- check_fun(superiority_reached, "sup", design_type)
    res_list[[which_cohort]]$Meta$fut_interim <- check_fun(futility_reached, "fut", design_type)
  } else {
    res_list[[which_cohort]]$Meta$sup_final_list <- superiority_reached
    res_list[[which_cohort]]$Meta$fut_final_list <- futility_reached
    res_list[[which_cohort]]$Meta$sup_final <- check_fun(superiority_reached, "sup", design_type)
    res_list[[which_cohort]]$Meta$fut_final <- check_fun(futility_reached, "fut", design_type)
  }

  # Write decisions into list
  if (interim) {
    if (res_list[[which_cohort]]$Meta$fut_interim) {
      res_list[[which_cohort]]$Meta$decision <- rep("STOP_FUT", 2)
    } else {
      if (res_list[[which_cohort]]$Meta$sup_interim) {
        res_list[[which_cohort]]$Meta$decision[1] <- "GO_SUP"
      } else {
        res_list[[which_cohort]]$Meta$decision[1] <- "CONTINUE"
      }
    }
  } else {
    if (res_list[[which_cohort]]$Meta$fut_final) {
      res_list[[which_cohort]]$Meta$decision[2] <- "STOP_FUT"
    } else {
      if (res_list[[which_cohort]]$Meta$sup_final) {
        res_list[[which_cohort]]$Meta$decision[2] <- "GO_SUP"
      } else {
        res_list[[which_cohort]]$Meta$decision[2] <- "STOP_N"
      }
    }
  }

  return(res_list)

}

