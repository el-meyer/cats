#' Checks whether decision criteria are met and updates trial results accordingly.
#'
#' Given a res_list object, checks the supplied decision criteria and saves the results in the res_list file.
#'
#' @param res_list         List item containing individual cohort trial results so far in a format used by the
#'                         other functions in this package
#'
#' @param which_cohort     Current cohort that should be evaluated
#'
#' @param test_strat       Testing strategy used; 1 = Combo vs. both monos, 2 = 1 + Add-on Mono vs. Placebo, 3 = 2 + Backbone mono vs. placebo
#'
#' @param sharing_type     What backbone and placebo data should be used for comparisons; Default is "all". Other options are "concurrent" or "dynamic" or "cohort".
#'
#' @param w                If dynamic borrowing, what is the prior choice for w. Default is 0.5.
#'
#' @param beta_prior       Prior parameter for all Beta Distributions. Default is 0.5.
#'
#' @param Bayes_Sup        List of matrices with rows corresponding to number of multiple Bayesian posterior two-arm combination criteria for superiority
#'
#' @param Bayes_Fut        List of matrices with rows corresponding to number of multiple Bayesian posterior two-arm combination criteria for futility
#'
#' @param Bayes_SA_Sup     List of matrices with rows corresponding to number of multiple Bayesian posterior single-arm combination criteria for superiority
#'
#' @param Bayes_SA_Fut     List of matrices with rows corresponding to number of multiple Bayesian posterior single-arm combination criteria for futility
#'
#' @param P_Sup            List with sublists corresponding to number of multiple frequentist test-based combination criteria for superiority
#'
#' @param P_Fut            List with sublists corresponding to number of multiple frequentist test-based combination criteria for futility
#'
#' @param Est_Sup_Fut      List with sublists corresponding to number of multiple point estimate based combination criteria for superiority and futility
#'
#' @param CI_Sup_Fut       List with sublists corresponding to number of multiple confidence interval based combination criteria for superiority and futility
#'
#' @param interim          Is the analysis conducted an interim or a final analysis?
#'
#' @param ...              Further arguments inherited from upper layer functions
#'
#' @return List containing original res_list and results of decision rules
#'
#' @examples
#'
#' # Example 1
#'
#' res_list <- list(c(list(decision = rep("none", 2), alloc_ratio = c(1,1,1,1),
#'                    n_thresh = c(Inf, 210)),
#'            rep(list(list(rr = NULL, resp_bio = NULL, resp_hist = NULL, n = NULL)), 4)))
#'
#' names(res_list)[1] <- paste0("Cohort", 1)
#' names(res_list[[1]])[4:7] <- c("Comb", "Mono", "Back", "Plac")
#' res_list[[1]][[4]]$rr <- 0.2
#' res_list[[1]][[5]]$rr <- 0.15
#' res_list[[1]][[6]]$rr <- 0.15
#' res_list[[1]][[7]]$rr <- 0.10
#'
#' r141 <- rbinom(1, 70, prob = res_list[[1]][[4]]$rr)
#' res_list[[1]][[4]]$resp_bio <- gtools::permute(c(rep(1, r141), rep(0, 70 - r141)))
#' r151 <- rbinom(1, 70, prob = res_list[[1]][[5]]$rr)
#' res_list[[1]][[5]]$resp_bio <- gtools::permute(c(rep(1, r151), rep(0, 70 - r151)))
#' r161 <- rbinom(1, 70, prob = res_list[[1]][[6]]$rr)
#' res_list[[1]][[6]]$resp_bio <- gtools::permute(c(rep(1, r161), rep(0, 70 - r161)))
#' r171 <- rbinom(1, 70, prob = res_list[[1]][[7]]$rr)
#' res_list[[1]][[7]]$resp_bio <- gtools::permute(c(rep(1, r171), rep(0, 70 - r171)))
#' r142 <- rbinom(1, 70, prob = res_list[[1]][[4]]$rr)
#' res_list[[1]][[4]]$resp_hist <- gtools::permute(c(rep(1, r142), rep(0, 70 - r142)))
#' r152 <- rbinom(1, 70, prob = res_list[[1]][[5]]$rr)
#' res_list[[1]][[5]]$resp_hist <- gtools::permute(c(rep(1, r152), rep(0, 70 - r152)))
#' r162 <- rbinom(1, 70, prob = res_list[[1]][[6]]$rr)
#' res_list[[1]][[6]]$resp_hist <- gtools::permute(c(rep(1, r162), rep(0, 70 - r162)))
#' r172 <- rbinom(1, 70, prob = res_list[[1]][[7]]$rr)
#' res_list[[1]][[7]]$resp_hist <- gtools::permute(c(rep(1, r172), rep(0, 70 - r172)))
#'
#' res_list[[1]][[4]]$n <- rep(1, 70)
#' res_list[[1]][[5]]$n <- rep(1, 70)
#' res_list[[1]][[6]]$n <- rep(1, 70)
#' res_list[[1]][[7]]$n <- rep(1, 70)
#'
#' # Comparison Combo vs Mono
#' Bayes_Sup1 <- matrix(nrow = 3, ncol = 3)
#' Bayes_Sup1[1,] <- c(0.00, 0.95, 0.90)
#' Bayes_Sup1[2,] <- c(0.10, 0.80, 0.75)
#' Bayes_Sup1[3,] <- c(0.15, 0.50, 1.00)
#' # Comparison Combo vs Backbone
#' Bayes_Sup2 <- matrix(nrow = 3, ncol = 3)
#' Bayes_Sup2[1,] <- c(0.00, 0.95, 0.90)
#' Bayes_Sup2[2,] <- c(NA, NA, NA)
#' Bayes_Sup2[3,] <- c(NA, NA, NA)
#' # Comparison Mono vs Placebo
#' Bayes_Sup3 <- matrix(nrow = 3, ncol = 3)
#' Bayes_Sup3[1,] <- c(0.00, 0.95, 0.90)
#' Bayes_Sup3[2,] <- c(0.10, 0.80, 0.75)
#' Bayes_Sup3[3,] <- c(NA, NA, NA)
#' #' # Comparison Backbone vs Placebo
#' Bayes_Sup4 <- matrix(nrow = 3, ncol = 3)
#' Bayes_Sup4[1,] <- c(0.00, 0.95, 0.90)
#' Bayes_Sup4[2,] <- c(0.10, 0.80, 0.75)
#' Bayes_Sup4[3,] <- c(NA, NA, NA)
#' Bayes_Sup <- list(list(Bayes_Sup1, Bayes_Sup2, Bayes_Sup3, Bayes_Sup4),
#'                   list(Bayes_Sup1, Bayes_Sup2, Bayes_Sup3, Bayes_Sup4))
#'
#' sharing_type <- "all"
#' interim <- TRUE
#' which_cohort <- 1
#' missing_prob <- 0.5
#' seed_missing <- 100
#'
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
#' # Combo
#' Bayes_SA_Sup1 <- matrix(nrow = 1, ncol = 3)
#' Bayes_SA_Sup1[1,] <- c(0.20, 0.95, 0.90)
#' # Mono
#' Bayes_SA_Sup2 <- matrix(nrow = 1, ncol = 3)
#' Bayes_SA_Sup2[1,] <- c(0.15, 0.80, 0.75)
#' # Backbone
#' Bayes_SA_Sup3 <- matrix(nrow = 1, ncol = 3)
#' Bayes_SA_Sup3[1,] <- c(0.15, 0.80, 0.75)
#' # Placebo
#' Bayes_SA_Sup4 <- matrix(nrow = 1, ncol = 3)
#' Bayes_SA_Sup4[1,] <- c(0.15, 0.80, 0.75)
#'
#' Bayes_SA_Sup <- list(list(Bayes_SA_Sup1, Bayes_SA_Sup2, Bayes_SA_Sup3, Bayes_SA_Sup4),
#'                      list(Bayes_SA_Sup1, Bayes_SA_Sup2, Bayes_SA_Sup3, Bayes_SA_Sup4))
#'
#' ## Combo
#' Bayes_SA_Fut1 <- matrix(nrow = 1, ncol = 2)
#' Bayes_SA_Fut1[1,] <- c(0.20, 0.50)
#' # Mono
#' Bayes_SA_Fut2 <- matrix(nrow = 1, ncol = 2)
#' Bayes_SA_Fut2[1,] <- c(0.15, 0.50)
#' # Backbone
#' Bayes_SA_Fut3 <- matrix(nrow = 1, ncol = 2)
#' Bayes_SA_Fut3[1,] <- c(0.15, 0.50)
#' # Placebo
#' Bayes_SA_Fut4 <- matrix(nrow = 1, ncol = 2)
#' Bayes_SA_Fut4[1,] <- c(0.15, 0.50)
#'
#' Bayes_SA_Fut <- list(list(Bayes_SA_Fut1, Bayes_SA_Fut2, Bayes_SA_Fut3, Bayes_SA_Fut4),
#'                      list(Bayes_SA_Fut1, Bayes_SA_Fut2, Bayes_SA_Fut3, Bayes_SA_Fut4))
#'
#' # Comparison Combo vs Mono
#' P_Sup1 <- list(list(
#' testfun = function(x) stats::prop.test(x, alternative = "less", correct = FALSE),
#' p_sup = 0.025, p_prom = 0.10, p_adj = "B"))
#' # Comparison Combo vs Backbone
#' P_Sup2 <- list(list(
#' testfun = function(x) stats::prop.test(x, alternative = "less", correct = FALSE),
#' p_sup = 0.025, p_prom = 0.10, p_adj = "B"))
#' # Comparison Mono vs Placebo
#' P_Sup3 <- list(list(
#' testfun = function(x) stats::prop.test(x, alternative = "less", correct = FALSE),
#' p_sup = 0.050, p_prom = 0.10, p_adj = "B"))
#' P_Sup4 <- list(list(
#' testfun = function(x) stats::prop.test(x, alternative = "less", correct = FALSE),
#' p_sup = 0.050, p_prom = 0.10, p_adj = "B"))
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
#' # Comparison Combo vs Mono
#' Est_Sup_Fut1 <- list(list(est = "AR", p_hat_sup = 0.6, p_hat_fut = 0.1, p_hat_prom = 0.5))
#' # Comparison Combo vs Backbone
#' Est_Sup_Fut2 <- list(list(est = "RR", p_hat_sup = 1.25, p_hat_fut = 0.75, p_hat_prom = 1.5))
#' # Comparison Mono vs Placebo
#' Est_Sup_Fut3 <- list(list(est = "OR", p_hat_sup = 1.50, p_hat_fut = 0.75, p_hat_prom = 2))
#' Est_Sup_Fut4 <- list(list(est = "OR", p_hat_sup = 1.50, p_hat_fut = 0.75, p_hat_prom = 2))
#' Est_Sup_Fut <- list(list(Est_Sup_Fut1, Est_Sup_Fut2, Est_Sup_Fut3, Est_Sup_Fut4),
#'                     list(Est_Sup_Fut1, Est_Sup_Fut2, Est_Sup_Fut3, Est_Sup_Fut4))
#'
#' # Comparison Combo vs Mono
#' CI_Sup_Fut1 <- list(list(est = "AR", ci = 0.95, p_hat_lower_sup = 0.35,
#'                    p_hat_upper_fut = 0.25, p_hat_lower_prom = 0.3))
#' # Comparison Combo vs Backbone
#' CI_Sup_Fut2 <- list(list(est = "RR", ci = 0.95, p_hat_lower_sup = 1.10,
#'                    p_hat_upper_fut = 1.10, p_hat_lower_prom = 1.05))
#' # Comparison Mono vs Placebo
#' CI_Sup_Fut3 <- list(list(est = "OR", ci = 0.95, p_hat_lower_sup = 1.20,
#'                    p_hat_upper_fut = 1.20, p_hat_lower_prom = 1.10))
#' CI_Sup_Fut4 <- list(list(est = "OR", ci = 0.95, p_hat_lower_sup = 1.20,
#'                    p_hat_upper_fut = 1.20, p_hat_lower_prom = 1.10))
#' CI_Sup_Fut <- list(list(CI_Sup_Fut1, CI_Sup_Fut2, CI_Sup_Fut3, CI_Sup_Fut4),
#'                    list(CI_Sup_Fut1, CI_Sup_Fut2, CI_Sup_Fut3, CI_Sup_Fut4))
#'
#' make_decision_trial(res_list = res_list, which_cohort = which_cohort, interim = interim,
#' Bayes_Sup = Bayes_Sup, sharing_type = sharing_type,
#' Bayes_Fut = Bayes_Fut, Bayes_SA_Sup = Bayes_SA_Sup, Bayes_SA_Fut = Bayes_SA_Fut, P_Sup = P_Sup,
#' P_Fut = P_Fut, Est_Sup_Fut = Est_Sup_Fut, CI_Sup_Fut = CI_Sup_Fut
#' )
#'
#' @export
make_decision_trial <- function(res_list, which_cohort, test_strat = 3, sharing_type = "all",
                                  Bayes_Sup = NULL, Bayes_Fut = NULL, Bayes_SA_Sup = NULL, Bayes_SA_Fut = NULL,
                                  w = 0.5, P_Sup = NULL, P_Fut = NULL, Est_Sup_Fut = NULL, CI_Sup_Fut = NULL,
                                  interim, beta_prior = 0.5, missing_prob = 0, seed_missing = 1, ...) {

  # todo: allow Estimates and CIs to be single arm rules for combo, mono and back
  # (atm only implemented as two-arm decision rules)

  ########## Parameter Explanation and Variable and Function Initiation ##########

  # Bayes_Sup = Matrix with rows corresponding to number of multiple Bayesian two-arm posterior combination criteria for superiority
  # First element is always the delta, second element is the GO-threshold and third is the promising-threshold
  # The rule is: If P(P_e > P_c + delta) > GO-threshold, then GO
  # The second rule is: If promising_threshold < P(P_e > P_c + delta) < GO-threshold, then declare promising

  # Bayes_Fut = Matrix with rows corresponding to number of multiple Bayesian two-arm posterior combination criteria for futility
  # First element is always the delta and the second element is the STOP-threshold
  # The rule is: If P(P_e > P_c + delta) < STOP-threshold, then STOP

  # Bayes_SA_Sup = Matrix with rows corresponding to number of multiple Bayesian single-arm posterior combination criteria for superiority
  # First element is always the margin, second element is the GO-threshold and third is the promising-threshold
  # The rule is: If P(P_e > margin) > GO-threshold, then GO
  # The second rule is: If promising_threshold < P(P_e > margin) < GO-threshold, then declare promising

  # Bayes_SA_Fut = Matrix with rows corresponding to number of multiple Bayesian single-arm posterior combination criteria for futility
  # First element is always the margin and the second element is the STOP-threshold
  # The rule is: If P(P_e > margin) < STOP-threshold, then STOP

  # P_Sup = List with sublists corresponding to number of multiple frequentist test-based combination criteria for superiority
  # Every sublist has the same structure:
  # The first element is the Test Function. This needs to be a function out of a list of allowed functions, or any
  # custom function that takes a 2x2 contingency table as input
  # The second element is the Superiority Threshold (i.e. GO, if p < p_sup)
  # The third element is the Promising Threshold (i.e. PROMISING, if p_sup < p < p_prom)
  # The fourth element contains the p-value multiplicity adjustment method at every decision

  # P_Fut = List with sublists corresponding to number of multiple frequentist test-based combination criteria for futility
  # Every sublist has the same structure:
  # The first element is the Test Function. This needs to be a function out of a list of allowed functions, or any
  # custom function that takes a 2x2 contingency table as input
  # The second element is the Superiority Threshold (i.e. STOP, if p > p_fut)

  # Est_Sup_Fut = List containing in sublists the rules corresponding to superiority and futility decision rules based on the point estimates
  # Every sublist has the same structure:
  # The first element is the type of point estimate. Choices are "AR", "RR", "OR".
  # AR estimation is naive
  # RR estimation via small sample adjustment UMLE from epitools package
  # OR using median-unbiased estimate from epitools package
  # The second element is the threshold required to declare superiority, e.g. GO if p_hat > p_hat_sup
  # The third element is the threshold required to declare futility, e.g. STOP if p_hat < p_hat_fut
  # The fourth element is the threshold required to declare promising, e.g. PROMISING if p_hat > p_hat_prom

  # CI_Sup_Fut = List containing in sublists the rules corresponding to superiority and futility decision rules based on confidence intervals
  # Every sublist has the same structure:
  # The first element is the type of confidence interval. Choices are "AR", "RR", "OR".
  # AR CI estimation via Clopper-Pearson
  # RR estimation via normal Approximation with small sample adjustment (Wald) epitools package
  # OR using mid-p exact CI from epitools package
  # The second element is the coverage probability of the confidence interval
  # The third element is the threshold required for the lower CI bound to declare superiority, e.g. GO if p_hat_lower > p_hat_lower_sup
  # The fourth element is the threshold required for the upper CI bound to declare futility, e.g. STOP if p_hat_upper < p_hat_upper_fut
  # The fifth element is the threshold required for the lower CI bound to declare promising, e.g. PROMISING if p_hat_lower > p_hat_lower_prom


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
  promising <- list()

  sup_reached <- fut_reached <- prom_reached <- NULL

  # Get responders and sample sizes
  resp_bio_tot <- resp_hist_tot <- n_tot <- NULL

  # Use temp dataset for number of responses and sample sizes (to allow for missing values)

  res_temp <- res_list

  # Combo
  if (interim) {
    resp_bio_tot[1] <- sum(res_temp[[which_cohort]]$Comb$resp_bio, na.rm = TRUE)
    n_tot[1] <- sum(res_temp[[which_cohort]]$Comb$n, na.rm = TRUE)
  } else {
    resp_hist_tot[1] <- sum(res_temp[[which_cohort]]$Comb$resp_hist_incl_missing, na.rm = TRUE)
    n_tot[1] <- sum(res_temp[[which_cohort]]$Comb$n_incl_missing, na.rm = TRUE)
  }

  # Mono
  if (interim) {
    resp_bio_tot[2] <- sum(res_temp[[which_cohort]]$Mono$resp_bio, na.rm = TRUE)
    n_tot[2] <- sum(res_temp[[which_cohort]]$Mono$n, na.rm = TRUE)
  } else {
    resp_hist_tot[2] <- sum(res_temp[[which_cohort]]$Mono$resp_hist_incl_missing, na.rm = TRUE)
    n_tot[2] <- sum(res_temp[[which_cohort]]$Mono$n_incl_missing, na.rm = TRUE)
  }

  # Plac
  if (length(res_temp[[which_cohort]]$alloc_ratio) == 4) {

    if (sharing_type == "cohort") {
      if (interim) {
        resp_bio_tot[4] <- sum(res_temp[[which_cohort]]$Plac$resp_bio, na.rm = TRUE)
        n_tot[4] <- sum(res_temp[[which_cohort]]$Plac$n, na.rm = TRUE)
      } else {
        resp_hist_tot[4] <- sum(res_temp[[which_cohort]]$Plac$resp_hist_incl_missing, na.rm = TRUE)
        n_tot[4] <- sum(res_temp[[which_cohort]]$Plac$n_incl_missing, na.rm = TRUE)
      }
    }

    if (sharing_type == "all") {
      if (interim) {
        resp_bio_tot[4] <- sum(sapply(res_temp, function(x) sum(x$Plac$resp_bio, na.rm = TRUE)), na.rm = TRUE)
        n_tot[4] <- sum(sapply(res_temp, function(x) sum(x$Plac$n, na.rm = TRUE)), na.rm = TRUE)
      } else {
        resp_hist_tot[4] <- sum(sapply(res_temp, function(x) sum(x$Plac$resp_hist_incl_missing, na.rm = TRUE)), na.rm = TRUE)
        n_tot[4] <- sum(sapply(res_temp, function(x) sum(x$Plac$n_incl_missing, na.rm = TRUE)), na.rm = TRUE)
      }
    }

    if (sharing_type == "concurrent") {
      conc <- which(!is.na(res_temp[[which_cohort]]$Plac$n))
      if (interim) {
        resp_bio_tot[4] <- sum(sapply(res_temp, function(x) sum(x$Plac$resp_bio[conc], na.rm = TRUE)), na.rm = TRUE)
        n_tot[4] <- sum(sapply(res_temp, function(x) sum(x$Plac$n[conc], na.rm = TRUE)), na.rm = TRUE)
      } else {
        resp_hist_tot[4] <- sum(sapply(res_temp, function(x) sum(x$Plac$resp_hist_incl_missing[conc], na.rm = TRUE)), na.rm = TRUE)
        n_tot[4] <- sum(sapply(res_temp, function(x) sum(x$Plac$n_incl_missing[conc], na.rm = TRUE)), na.rm = TRUE)
      }
    }

    if (sharing_type == "dynamic") {
      # if all are concurrent, no borrowing needed
      if (all(!is.na(res_temp[[which_cohort]]$Plac$n))) {
        if (interim) {
          resp_bio_tot[4] <- sum(sapply(res_temp, function(x) sum(x$Plac$resp_bio, na.rm = TRUE)), na.rm = TRUE)
          n_tot[4] <- sum(sapply(res_temp, function(x) sum(x$Plac$n, na.rm = TRUE)), na.rm = TRUE)
        } else {
          resp_hist_tot[4] <- sum(sapply(res_temp, function(x) sum(x$Plac$resp_hist_incl_missing, na.rm = TRUE)), na.rm = TRUE)
          n_tot[4] <- sum(sapply(res_temp, function(x) sum(x$Plac$n_incl_missing, na.rm = TRUE)), na.rm = TRUE)
        }
      } else {
        # do borrowing

        if (interim) {
          # compute responders from cohort
          suc_bio_c <- sum(res_temp[[which_cohort]]$Plac$resp_bio, na.rm = TRUE)
          N_c <- sum(res_temp[[which_cohort]]$Plac$n, na.rm = TRUE)

          # compute historical responders
          suc_bio_h <- sum(sapply(res_temp, function(x) sum(x$Plac$resp_bio, na.rm = TRUE)), na.rm = TRUE) - suc_bio_c
          N_h <- sum(sapply(res_temp, function(x) sum(x$Plac$n, na.rm = TRUE)), na.rm = TRUE) - N_c

          w1_bio <- w * beta(suc_bio_c + suc_bio_h + beta_prior, N_h + N_c - suc_bio_c - suc_bio_h + beta_prior) /
            beta(suc_bio_h + beta_prior, N_h - suc_bio_h + beta_prior)
          w2_bio <- (1 - w) * beta(suc_bio_c + beta_prior, N_c - suc_bio_c + beta_prior) / beta(beta_prior, beta_prior)
          w1n_bio <- w1_bio / (w1_bio + w2_bio)
          w2n_bio <- w2_bio / (w1_bio + w2_bio)

          m_a_bio <- w1n_bio * (suc_bio_c + suc_bio_h + beta_prior) + w2n_bio * (suc_bio_c + beta_prior)
          m_b_bio <- w1n_bio * (N_h + N_c - suc_bio_c - suc_bio_h + beta_prior) + w2n_bio * (N_c - suc_bio_c + beta_prior)

          n_tot[4] <- round(m_a_bio + m_b_bio)
          resp_bio_tot[4] <- round(m_a_bio)

        } else {

          suc_hist_c <- sum(res_temp[[which_cohort]]$Plac$resp_hist_incl_missing, na.rm = TRUE)
          N_c <- sum(res_temp[[which_cohort]]$Plac$n_incl_missing, na.rm = TRUE)

          # compute historical responders
          suc_hist_h <- sum(sapply(res_temp, function(x) sum(x$Plac$resp_hist_incl_missing, na.rm = TRUE)), na.rm = TRUE) - suc_hist_c
          N_h <- sum(sapply(res_temp, function(x) sum(x$Plac$n_incl_missing, na.rm = TRUE)), na.rm = TRUE) - N_c

          w1_hist <- w * beta(suc_hist_c + suc_hist_h + beta_prior, N_h + N_c - suc_hist_c - suc_hist_h + beta_prior) /
            beta(suc_hist_h + beta_prior, N_h - suc_hist_h + beta_prior)
          w2_hist <- (1 - w) * beta(suc_hist_c + beta_prior, N_c - suc_hist_c + beta_prior) / beta(beta_prior, beta_prior)
          w1n_hist <- w1_hist / (w1_hist + w2_hist)
          w2n_hist <- w2_hist / (w1_hist + w2_hist)

          m_a_hist <- w1n_hist * (suc_hist_c + suc_hist_h + beta_prior) + w2n_hist * (suc_hist_c + beta_prior)
          m_b_hist <- w1n_hist * (N_h + N_c - suc_hist_c - suc_hist_h + beta_prior) + w2n_hist * (N_c - suc_hist_c + beta_prior)

          n_tot[4] <- round(m_a_hist + m_b_hist)
          resp_hist_tot[4] <- round(m_a_hist)

        }
      }
    }
  }

  # separately for backbone

  if (sharing_type == "cohort") {
    if (interim) {
      resp_bio_tot[3] <- sum(res_temp[[which_cohort]]$Back$resp_bio, na.rm = TRUE)
      n_tot[3] <- sum(res_temp[[which_cohort]]$Back$n, na.rm = TRUE)
    } else {
      resp_hist_tot[3] <- sum(res_temp[[which_cohort]]$Back$resp_hist_incl_missing, na.rm = TRUE)
      n_tot[3] <- sum(res_temp[[which_cohort]]$Back$n_incl_missing, na.rm = TRUE)
    }

  }

  if (sharing_type == "all") {
    if (interim) {
      resp_bio_tot[3] <- sum(sapply(res_temp, function(x) sum(x$Back$resp_bio, na.rm = TRUE)), na.rm = TRUE)
      n_tot[3] <- sum(sapply(res_temp, function(x) sum(x$Back$n, na.rm = TRUE)), na.rm = TRUE)
    } else {
      resp_hist_tot[3] <- sum(sapply(res_temp, function(x) sum(x$Back$resp_hist_incl_missing, na.rm = TRUE)), na.rm = TRUE)
      n_tot[3] <- sum(sapply(res_temp, function(x) sum(x$Back$n_incl_missing, na.rm = TRUE)), na.rm = TRUE)
    }
  }

  if (sharing_type == "concurrent") {
    conc <- which(!is.na(res_temp[[which_cohort]]$Back$n))
    if (interim) {
      resp_bio_tot[3] <- sum(sapply(res_temp, function(x) sum(x$Back$resp_bio[conc], na.rm = TRUE)), na.rm = TRUE)
      n_tot[3] <- sum(sapply(res_temp, function(x) sum(x$Back$n[conc], na.rm = TRUE)), na.rm = TRUE)
    } else {
      resp_hist_tot[3] <- sum(sapply(res_temp, function(x) sum(x$Back$resp_hist_incl_missing[conc], na.rm = TRUE)), na.rm = TRUE)
      n_tot[3] <- sum(sapply(res_temp, function(x) sum(x$Back$n_incl_missing[conc], na.rm = TRUE)), na.rm = TRUE)
    }

  }

  if (sharing_type == "dynamic") {
    # if all are concurrent, no borrowing needed
    if (all(!is.na(res_temp[[which_cohort]]$Back$n))) {
      if (interim) {
        resp_bio_tot[3] <- sum(sapply(res_temp, function(x) sum(x$Back$resp_bio, na.rm = TRUE)), na.rm = TRUE)
        n_tot[3] <- sum(sapply(res_temp, function(x) sum(x$Back$n, na.rm = TRUE)), na.rm = TRUE)
      } else {
        resp_hist_tot[3] <- sum(sapply(res_temp, function(x) sum(x$Back$resp_hist_incl_missing, na.rm = TRUE)), na.rm = TRUE)
        n_tot[3] <- sum(sapply(res_temp, function(x) sum(x$Back$n_incl_missing, na.rm = TRUE)), na.rm = TRUE)
      }

    } else {
      # do borrowing

      if (interim) {
        # compute responders from cohort
        suc_bio_c <- sum(res_temp[[which_cohort]]$Back$resp_bio, na.rm = TRUE)
        N_c <- sum(res_temp[[which_cohort]]$Back$n, na.rm = TRUE)

        # compute historical responders
        suc_bio_h <- sum(sapply(res_temp, function(x) sum(x$Back$resp_bio, na.rm = TRUE)), na.rm = TRUE) - suc_bio_c
        N_h <- sum(sapply(res_temp, function(x) sum(x$Back$n, na.rm = TRUE)), na.rm = TRUE) - N_c

        w1_bio <- w * beta(suc_bio_c + suc_bio_h + beta_prior, N_h + N_c - suc_bio_c - suc_bio_h + beta_prior) /
          beta(suc_bio_h + beta_prior, N_h - suc_bio_h + beta_prior)
        w2_bio <- (1 - w) * beta(suc_bio_c + beta_prior, N_c - suc_bio_c + beta_prior) / beta(beta_prior, beta_prior)
        w1n_bio <- w1_bio / (w1_bio + w2_bio)
        w2n_bio <- w2_bio / (w1_bio + w2_bio)

        m_a_bio <- w1n_bio * (suc_bio_c + suc_bio_h + beta_prior) + w2n_bio * (suc_bio_c + beta_prior)
        m_b_bio <- w1n_bio * (N_h + N_c - suc_bio_c - suc_bio_h + beta_prior) + w2n_bio * (N_c - suc_bio_c + beta_prior)

        n_tot[3] <- round(m_a_bio + m_b_bio)
        resp_bio_tot[3] <- round(m_a_bio)

      } else {

        # compute responders from cohort
        suc_hist_c <- sum(res_temp[[which_cohort]]$Back$resp_hist_incl_missing, na.rm = TRUE)
        N_c <- sum(res_temp[[which_cohort]]$Back$n_incl_missing, na.rm = TRUE)

        # compute historical responders
        suc_hist_h <- sum(sapply(res_temp, function(x) sum(x$Back$resp_hist_incl_missing, na.rm = TRUE)), na.rm = TRUE) - suc_hist_c
        N_h <- sum(sapply(res_temp, function(x) sum(x$Back$n_incl_missing, na.rm = TRUE)), na.rm = TRUE) - N_c

        w1_hist <- w * beta(suc_hist_c + suc_hist_h + beta_prior, N_h + N_c - suc_hist_c - suc_hist_h + beta_prior) /
          beta(suc_hist_h + beta_prior, N_h - suc_hist_h + beta_prior)
        w2_hist <- (1 - w) * beta(suc_hist_c + beta_prior, N_c - suc_hist_c + beta_prior) / beta(beta_prior, beta_prior)
        w1n_hist <- w1_hist / (w1_hist + w2_hist)
        w2n_hist <- w2_hist / (w1_hist + w2_hist)

        m_a_hist <- w1n_hist * (suc_hist_c + suc_hist_h + beta_prior) + w2n_hist * (suc_hist_c + beta_prior)
        m_b_hist <- w1n_hist * (N_h + N_c - suc_hist_c - suc_hist_h + beta_prior) + w2n_hist * (N_c - suc_hist_c + beta_prior)

        n_tot[3] <- round(m_a_hist + m_b_hist)
        resp_hist_tot[3] <- round(m_a_hist)
      }


    }
  }

  if (interim) {
    resp_tot <- resp_bio_tot
  } else {
    resp_tot <- resp_hist_tot
  }

  j_fin <- test_strat + 1

  if (length(res_list[[which_cohort]]$alloc_ratio) == 3) {
    j_fin <- 2
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
      assign(paste0("prob_sup_vec", k), rep(NA, 4))
    }

    assign(paste0("sup_reached"), matrix(NA, nrow = max(sapply(Bayes_Sup, function(x) nrow(x))), ncol = 4))
    assign(paste0("prom_reached"), matrix(NA, nrow = max(sapply(Bayes_Sup, function(x) nrow(x))), ncol = 4))

    # For each of the comparisons separately, evaluate the applying decision rules
    # j=1: Combo vs Mono, j=2: Combo vs Back, j=3: Mono vs Placebo

    for (i in 1:max(sapply(Bayes_Sup, function(x) nrow(x)))) {
      for (j in 1:j_fin) {
        # Get which type of comparison should be conducted
        if (j == 1) {
          c1 <- 1
          c2 <- 2
        }
        if (j == 2) {
          c1 <- 1
          c2 <- 3
        }
        if (j == 3) {
          c1 <- 2
          c2 <- 4
        }
        if (j == 4) {
          c1 <- 3
          c2 <- 4
        }

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

        # Check whether promising reached
        eval(parse(text = paste(paste0("prom_reached[", i, ",", j, "]"), "<-",
                                get(paste0("prob_sup_vec", i))[j] <= Bayes_Sup[[j]][i,2] &
                                  get(paste0("prob_sup_vec", i))[j] >= Bayes_Sup[[j]][i,3])))
      }
      # Write results in res_list
      eval(parse(text = paste(paste0("res_list[[which_cohort]]$prob_sup", i, "<-",
                                     "c(res_list[[which_cohort]]$prob_sup", i, ",prob_sup_vec", i, ")"))))
    }

    # Add decisions from all decision rules for each arm together
    superiority_reached <- c(superiority_reached, list(sup_reached))
    promising <- c(promising, list(prom_reached))

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
      assign(paste0("prob_fut_vec", k), rep(NA, 4))
    }

    assign(paste0("fut_reached"), matrix(NA, nrow = max(sapply(Bayes_Fut, function(x) nrow(x))), ncol = 4))

    # For each of the comparisons separately, evaluate the applying decision rules
    # j=1: Combo vs Mono, j=2: Combo vs Back, j=3: Mono vs Placebo

    for (i in 1:max(sapply(Bayes_Fut, function(x) nrow(x)))) {
      for (j in 1:j_fin) {
        # Get which type of comparison should be conducted
        if (j == 1) {
          c1 <- 1
          c2 <- 2
        }
        if (j == 2) {
          c1 <- 1
          c2 <- 3
        }
        if (j == 3) {
          c1 <- 2
          c2 <- 4
        }
        if (j == 4) {
          c1 <- 3
          c2 <- 4
        }
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
      eval(parse(text = paste(paste0("res_list[[which_cohort]]$prob_fut", i, "<-",
                                     "c(res_list[[which_cohort]]$prob_fut", i, ",prob_fut_vec", i, ")"))))
    }

    # Add decisions from all decision rules for each arm together
    futility_reached <- c(futility_reached, list(fut_reached))

  }

  ########## Bayesian Single-Arm Superiority Criteria ###############

  if (!is.null(Bayes_SA_Sup)) {
    # Get the decision rule for interim/final
    if (interim) {
      Bayes_SA_Sup <- Bayes_SA_Sup[[1]]
    } else {
      Bayes_SA_Sup <- Bayes_SA_Sup[[2]]
    }

    # Firstly initiate as many vectors as will maximum be needed (which is determined by the most complex decision rule)
    for (k in 1:max(sapply(Bayes_SA_Sup, function(x) nrow(x)))) {
      assign(paste0("prob_sup_sa_vec", k), rep(NA, 4))
    }
    assign(paste0("sup_reached"), matrix(NA, nrow = max(sapply(Bayes_SA_Sup, function(x) nrow(x))), ncol = 4))
    assign(paste0("prom_reached"), matrix(NA, nrow = max(sapply(Bayes_SA_Sup, function(x) nrow(x))), ncol = 4))

    # For each of the comparisons separately, evaluate the applying decision rules
    # j=1: Combo vs Mono, j=2: Combo vs Back, j=3: Mono vs Placebo

    for (i in 1:max(sapply(Bayes_SA_Sup, function(x) nrow(x)))) {
      for (j in 1:j_fin) {
        # Get which type of comparison should be conducted
        # Since this is SA, the first rules apply to combo, the second rules to mono and the third to backbone
        if (j == 1) {
          c <- 1
        }
        if (j == 2) {
          c <- 2
        }
        if (j == 3) {
          c <- 3
        }
        if (j == 4) {
          c <- 4
        }

        # Evaluate decision rule, but if non-existant, give NA
        if (is.na(Bayes_SA_Sup[[j]][i,1])) {
          eval(parse(text = paste(paste0("prob_sup_sa_vec", i, "[", j, "]"), "<- NA")))
        } else {
          eval(parse(text = paste(paste0("prob_sup_sa_vec", i, "[", j, "]"), "<-",
                                  1 - stats::pbeta(Bayes_SA_Sup[[j]][i,1],
                                            beta_prior + resp_tot[c],
                                            beta_prior + n_tot[c] - resp_tot[c]))))
        }

        # Check whether superiority reached
        eval(parse(text = paste(paste0("sup_reached[", i, ",", j, "]"), "<-",
                                get(paste0("prob_sup_sa_vec", i))[j] > Bayes_SA_Sup[[j]][i,2])))

        # Check whether promising reached
        eval(parse(text = paste(paste0("prom_reached[", i, ",", j, "]"), "<-",
                                get(paste0("prob_sup_sa_vec", i))[j] <= Bayes_SA_Sup[[j]][i,2] &
                                  get(paste0("prob_sup_sa_vec", i))[j] >= Bayes_SA_Sup[[j]][i,3])))
      }
      # Write results in res_list
      eval(parse(text = paste(paste0("res_list[[which_cohort]]$prob_sup_sa", i, "<-",
                                     "c(res_list[[which_cohort]]$prob_sup_sa", i, ",prob_sup_sa_vec", i, ")"))))
    }

    # Add decisions from all decision rules for each arm together
    superiority_reached <- c(superiority_reached, list(sup_reached))
    promising <- c(promising, list(prom_reached))

  }


  ########## Bayesian Single-Arm Futility Criteria ##########

  if (!is.null(Bayes_SA_Fut)) {
    # Get the decision rule for interim/final
    if (interim) {
      Bayes_SA_Fut <- Bayes_SA_Fut[[1]]
    } else {
      Bayes_SA_Fut <- Bayes_SA_Fut[[2]]
    }

    # Firstly initiate as many vectors as will maximum be needed (which is determined by the most complex decision rule)
    for (k in 1:max(sapply(Bayes_SA_Fut, function(x) nrow(x)))) {
      assign(paste0("prob_fut_sa_vec", k), rep(NA, 4))
    }
    assign(paste0("fut_reached"), matrix(NA, nrow = max(sapply(Bayes_SA_Fut, function(x) nrow(x))), ncol = 4))

    # For each of the comparisons separately, evaluate the applying decision rules
    # j=1: Combo vs Mono, j=2: Combo vs Back, j=3: Mono vs Placebo

    for (i in 1:max(sapply(Bayes_SA_Fut, function(x) nrow(x)))) {
      for (j in 1:j_fin) {
        # Get which type of comparison should be conducted
        # Since this is SA, the first rules apply to combo, the second rules to mono and the third to backbone
        if (j == 1) {
          c <- 1
        }
        if (j == 2) {
          c <- 2
        }
        if (j == 3) {
          c <- 3
        }
        if (j == 4) {
          c <- 4
        }

        # Evaluate decision rule, but if non-existant, give NA
        if (is.na(Bayes_SA_Fut[[j]][i,1])) {
          eval(parse(text = paste(paste0("prob_fut_sa_vec", i, "[", j, "]"), "<- NA")))
        } else {
          eval(parse(text = paste(paste0("prob_fut_sa_vec", i, "[", j, "]"), "<-",
                                  1 - stats::pbeta(Bayes_SA_Fut[[j]][i,1],
                                            beta_prior + resp_tot[c],
                                            beta_prior + n_tot[c] - resp_tot[c]))))
        }

        # Check whether futility reached
        eval(parse(text = paste(paste0("fut_reached[", i, ",", j, "]"), "<-",
                                get(paste0("prob_fut_sa_vec", i))[j] < Bayes_SA_Fut[[j]][i,2])))

      }
      # Write results in res_list
      eval(parse(text = paste(paste0("res_list[[which_cohort]]$prob_fut_sa", i, "<-",
                                     "c(res_list[[which_cohort]]$prob_fut_sa", i, ",prob_fut_sa_vec", i, ")"))))
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
      assign(paste0("p_values", k), rep(NA, 4))
      assign(paste0("p_adjusted", k), rep(NA, 4))
    }

    assign(paste0("sup_reached"), matrix(NA, nrow = max(sapply(P_Sup, function(x) length(x))), ncol = 4))
    assign(paste0("prom_reached"), matrix(NA, nrow = max(sapply(P_Sup, function(x) length(x))), ncol = 4))

    # For each of the comparisons separately, evaluate the applying decision rules
    # j=1: Combo vs Mono, j=2: Combo vs Back, j=3: Mono vs Placebo

    for (i in 1:max(sapply(P_Sup, function(x) length(x)))) {
      for (j in 1:j_fin) {
        # Get which type of comparison should be conducted
        if (j == 1) {
          c1 <- 1
          c2 <- 2
        }
        if (j == 2) {
          c1 <- 1
          c2 <- 3
        }
        if (j == 3) {
          c1 <- 2
          c2 <- 4
        }
        if (j == 4) {
          c1 <- 3
          c2 <- 4
        }

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

        # Check whether promising reached
        eval(parse(text = paste(paste0("prom_reached[", i, ",", j, "]"), "<-",
                                get(paste0("p_adjusted", i))[j] <= P_Sup[[j]][[i]]$p_prom)))
      }
      # Write results in res_list
      eval(parse(text = paste(paste0("res_list[[which_cohort]]$p_values_sup", i, "<-",
                                     "c(res_list[[which_cohort]]$p_values_sup", i, ",p_adjusted", i, ")"))))
    }

    # Add decisions from all decision rules for each arm together
    superiority_reached <- c(superiority_reached, list(sup_reached))
    promising <- c(promising, list(prom_reached))

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
      assign(paste0("p_values", k), rep(NA, 4))
      assign(paste0("p_adjusted", k), rep(NA, 4))
    }

    assign(paste0("fut_reached"), matrix(NA, nrow = max(sapply(P_Fut, function(x) length(x))), ncol = 4))

    # For each of the comparisons separately, evaluate the applying decision rules
    # j=1: Combo vs Mono, j=2: Combo vs Back, j=3: Mono vs Placebo

    for (i in 1:max(sapply(P_Fut, function(x) length(x)))) {
      for (j in 1:j_fin) {
        # Get which type of comparison should be conducted
        if (j == 1) {
          c1 <- 1
          c2 <- 2
        }
        if (j == 2) {
          c1 <- 1
          c2 <- 3
        }
        if (j == 3) {
          c1 <- 2
          c2 <- 4
        }
        if (j == 4) {
          c1 <- 3
          c2 <- 4
        }

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
      eval(parse(text = paste(paste0("res_list[[which_cohort]]$p_values_fut", i, "<-",
                                     "c(res_list[[which_cohort]]$p_values_fut", i, ",p_adjusted", i, ")"))))
    }

    # Add decisions from all decision rules for each arm together
    futility_reached <- c(futility_reached, list(fut_reached))

  }

  ########## Point Estimate Superiority/Futility Criteria ##########

  if (!is.null(Est_Sup_Fut)) {
    # Get the decision rule for interim/final
    if (interim) {
      Est_Sup_Fut <- Est_Sup_Fut[[1]]
    } else {
      Est_Sup_Fut <- Est_Sup_Fut[[2]]
    }

    # Firstly initiate as many vectors as will maximum be needed (which is determined by the most complex decision rule)
    for (k in 1:max(sapply(Est_Sup_Fut, function(x) length(x)))) {
      assign(paste0("p_hat", k), rep(NA, 4))
    }

    assign(paste0("sup_reached"), matrix(NA, nrow = max(sapply(Est_Sup_Fut, function(x) length(x))), ncol = 4))
    assign(paste0("fut_reached"), matrix(NA, nrow = max(sapply(Est_Sup_Fut, function(x) length(x))), ncol = 4))
    assign(paste0("prom_reached"), matrix(NA, nrow = max(sapply(Est_Sup_Fut, function(x) length(x))), ncol = 4))

    # For each of the comparisons separately, evaluate the applying decision rules
    # j=1: Combo vs Mono, j=2: Combo vs Back, j=3: Mono vs Placebo

    for (i in 1:max(sapply(Est_Sup_Fut, function(x) length(x)))) {
      for (j in 1:j_fin) {
        # Get which type of comparison should be conducted
        if (j == 1) {
          c1 <- 1
          c2 <- 2
        }
        if (j == 2) {
          c1 <- 1
          c2 <- 3
        }
        if (j == 3) {
          c1 <- 2
          c2 <- 4
        }
        if (j == 4) {
          c1 <- 3
          c2 <- 4
        }

        # Evaluate decision rule, but if non-existant, give NA
        if (is.na(Est_Sup_Fut[[j]][[i]]$est)) {
          eval(parse(text = paste(paste0("p_hat", i, "[", j, "]"), "<- NA")))
        } else {

          if (Est_Sup_Fut[[j]][[i]]$est == "AR") {
            eval(parse(text = paste(paste0("p_hat", i, "[", j, "]"), "<-", resp_tot[j] / n_tot[j])))

          } else {
            if (Est_Sup_Fut[[j]][[i]]$est == "RR") {
              f <- match.fun(function(x) epitools::riskratio.small(x, rev = "columns")$measure[2,1])
            }
            if (Est_Sup_Fut[[j]][[i]]$est == "OR") {
              f <- match.fun(function(x) epitools::oddsratio(x, rev = "columns")$measure[2,1])
            }

            crosstab <-  matrix(c(resp_tot[c2], n_tot[c2] - resp_tot[c2], resp_tot[c1], n_tot[c1] - resp_tot[c1]),
                                nrow = 2, ncol = 2, byrow = TRUE)
            eval(parse(text = paste(paste0("p_hat", i, "[", j, "]"), "<-", f(crosstab))))
          }

        }

        # Check whether superiority reached
        eval(parse(text = paste(paste0("sup_reached[", i, ",", j, "]"), "<-",
                                get(paste0("p_hat", i))[j] > Est_Sup_Fut[[j]][[i]]$p_hat_sup)))

        # Check whether futility reached
        eval(parse(text = paste(paste0("fut_reached[", i, ",", j, "]"), "<-",
                                get(paste0("p_hat", i))[j] < Est_Sup_Fut[[j]][[i]]$p_hat_fut)))

        # Check whether promising reached
        eval(parse(text = paste(paste0("prom_reached[", i, ",", j, "]"), "<-",
                                get(paste0("p_hat", i))[j] >= Est_Sup_Fut[[j]][[i]]$p_hat_prom)))
      }
      # Write results in res_list
      eval(parse(text = paste(paste0("res_list[[which_cohort]]$p_hat", i, "<-",
                                     "c(res_list[[which_cohort]]$p_hat", i, ",p_hat", i, ")"))))
    }

    # Add decisions from all decision rules for each arm together
    superiority_reached <- c(superiority_reached, list(sup_reached))
    futility_reached <- c(futility_reached, list(fut_reached))
    promising <- c(promising, list(prom_reached))

  }


  ########## Confidence Interval Superiority/Futility Criteria ##########

  if (!is.null(CI_Sup_Fut)) {
    # Get the decision rule for interim/final
    if (interim) {
      CI_Sup_Fut <- CI_Sup_Fut[[1]]
    } else {
      CI_Sup_Fut <- CI_Sup_Fut[[2]]
    }

    # Firstly initiate as many vectors as will maximum be needed (which is determined by the most complex decision rule)
    for (k in 1:max(sapply(CI_Sup_Fut, function(x) length(x)))) {
      assign(paste0("p_hat_upper", k), rep(NA, 4))
      assign(paste0("p_hat_lower", k), rep(NA, 4))
    }

    assign(paste0("sup_reached"), matrix(NA, nrow = max(sapply(CI_Sup_Fut, function(x) length(x))), ncol = 4))
    assign(paste0("fut_reached"), matrix(NA, nrow = max(sapply(CI_Sup_Fut, function(x) length(x))), ncol = 4))
    assign(paste0("prom_reached"), matrix(NA, nrow = max(sapply(CI_Sup_Fut, function(x) length(x))), ncol = 4))

    # For each of the comparisons separately, evaluate the applying decision rules
    # j=1: Combo vs Mono, j=2: Combo vs Back, j=3: Mono vs Placebo

    for (i in 1:max(sapply(CI_Sup_Fut, function(x) length(x)))) {
      for (j in 1:j_fin) {
        # Get which type of comparison should be conducted
        if (j == 1) {
          c1 <- 1
          c2 <- 2
        }
        if (j == 2) {
          c1 <- 1
          c2 <- 3
        }
        if (j == 3) {
          c1 <- 2
          c2 <- 4
        }
        if (j == 4) {
          c1 <- 3
          c2 <- 4
        }

        # Evaluate decision rule, but if non-existant, give NA
        if (is.na(CI_Sup_Fut[[j]][[i]]$est)) {
          eval(parse(text = paste(paste0("p_hat", i, "[", j, "]"), "<- NA")))
        } else {

          if (CI_Sup_Fut[[j]][[i]]$est == "AR") {
            eval(parse(text = paste(paste0("p_hat_upper", i, "[", j, "]"), "<-", stats::binom.test(c(resp_tot[j], n_tot[j]))$conf.int[2])))
            eval(parse(text = paste(paste0("p_hat_lower", i, "[", j, "]"), "<-", stats::binom.test(c(resp_tot[j], n_tot[j]))$conf.int[1])))

          } else {
            if (CI_Sup_Fut[[j]][[i]]$est == "RR") {
              f_up <- match.fun(function(x) epitools::riskratio.small(x, conf.level = CI_Sup_Fut[[j]][[i]]$ci, rev = "columns")$measure[2,3])
              f_lo <- match.fun(function(x) epitools::riskratio.small(x, conf.level = CI_Sup_Fut[[j]][[i]]$ci, rev = "columns")$measure[2,2])
            }
            if (CI_Sup_Fut[[j]][[i]]$est == "OR") {
              f_up <- match.fun(function(x) epitools::oddsratio(x, conf.level = CI_Sup_Fut[[j]][[i]]$ci, rev = "columns")$measure[2,3])
              f_lo <- match.fun(function(x) epitools::oddsratio(x, conf.level = CI_Sup_Fut[[j]][[i]]$ci, rev = "columns")$measure[2,2])
            }

            crosstab <-  matrix(c(resp_tot[c2], n_tot[c2] - resp_tot[c2], resp_tot[c1], n_tot[c1] - resp_tot[c1]),
                                nrow = 2, ncol = 2, byrow = TRUE)
            eval(parse(text = paste(paste0("p_hat_upper", i, "[", j, "]"), "<-", f_up(crosstab))))
            eval(parse(text = paste(paste0("p_hat_lower", i, "[", j, "]"), "<-", f_lo(crosstab))))
          }
        }

        # Check whether superiority reached
        eval(parse(text = paste(paste0("sup_reached[", i, ",", j, "]"), "<-",
                                get(paste0("p_hat_lower", i))[j] > CI_Sup_Fut[[j]][[i]]$p_hat_lower_sup)))

        # Check whether futility reached
        eval(parse(text = paste(paste0("fut_reached[", i, ",", j, "]"), "<-",
                                get(paste0("p_hat_upper", i))[j] < CI_Sup_Fut[[j]][[i]]$p_hat_upper_fut)))

        # Check whether promising reached
        eval(parse(text = paste(paste0("prom_reached[", i, ",", j, "]"), "<-",
                                get(paste0("p_hat_lower", i))[j] >= CI_Sup_Fut[[j]][[i]]$p_hat_lower_prom)))
      }
      # Write results in res_list
      eval(parse(text = paste(paste0("res_list[[which_cohort]]$p_hat_upper", i, "<-",
                                     "c(res_list[[which_cohort]]$p_hat_upper", i, ",p_hat_upper", i, ")"))))
      eval(parse(text = paste(paste0("res_list[[which_cohort]]$p_hat_lower", i, "<-",
                                     "c(res_list[[which_cohort]]$p_hat_lower", i, ",p_hat_lower", i, ")"))))
    }

    # Add decisions from all decision rules for each arm together
    superiority_reached <- c(superiority_reached, list(sup_reached))
    futility_reached <- c(futility_reached, list(fut_reached))
    promising <- c(promising, list(prom_reached))

  }


  ########## Combine Results And Return List ##########

  # Combine results: Choose superiority, futility and evaluating according to rules specified
  # Look at futility first (also random safety stop probability)! Then at superiority, then at promising.

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
  if (interim) {
    res_list[[which_cohort]]$sup_interim_list <- superiority_reached
    res_list[[which_cohort]]$fut_interim_list <- futility_reached
    res_list[[which_cohort]]$prom_interim_list <- promising
    res_list[[which_cohort]]$sup_interim <- check_fun(superiority_reached, "sup")
    res_list[[which_cohort]]$fut_interim <- check_fun(futility_reached, "fut")
    res_list[[which_cohort]]$prom_interim <- check_fun(promising, "prom")
  } else {
    res_list[[which_cohort]]$sup_final_list <- superiority_reached
    res_list[[which_cohort]]$fut_final_list <- futility_reached
    res_list[[which_cohort]]$prom_final_list <- promising
    res_list[[which_cohort]]$sup_final <- check_fun(superiority_reached, "sup")
    res_list[[which_cohort]]$fut_final <- check_fun(futility_reached, "fut")
    res_list[[which_cohort]]$prom_final <- check_fun(promising, "prom")
  }

  # Write decisions into list
  if (interim) {
    if (res_list[[which_cohort]]$fut_interim) {
      res_list[[which_cohort]]$decision <- rep("STOP_FUT", 2)
    } else {
      if (res_list[[which_cohort]]$sup_interim) {
        res_list[[which_cohort]]$decision[1] <- "GO_SUP"
      } else {
        if (res_list[[which_cohort]]$prom_interim) {
          res_list[[which_cohort]]$decision[1] <- "PROMISING"
        } else {
          res_list[[which_cohort]]$decision[1] <- "CONTINUE"
        }
      }
    }
  } else {
    if (res_list[[which_cohort]]$fut_final) {
      res_list[[which_cohort]]$decision[2] <- "STOP_FUT"
    } else {
      if (res_list[[which_cohort]]$sup_final) {
        res_list[[which_cohort]]$decision[2] <- "GO_SUP"
      } else {
        res_list[[which_cohort]]$decision[2] <- "STOP_N"
      }
    }
  }

  return(res_list)

}

