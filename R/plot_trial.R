#' Plots the cohort trial study overview given stage data.
#'
#' Given a res_list object, plots things like final study design, indicating which arms were discontinued after how many patients etc..
#'
#' @param res_list    List item containing trial results so far in a format used by the other functions in this package
#'
#' @param unit        What is unit of observation in response rate plots: N_cohort or N_total?
#'
#' @examples
#'
#' random <- TRUE
#'
#' rr_comb <- c(1)
#' prob_comb_rr <- c(1)
#' rr_mono <- c(1,2)
#' prob_mono_rr <- c(0.2, 0.8)
#' rr_back <- c(2)
#' prob_back_rr <- c(1)
#' rr_plac <- c(0.10)
#' prob_plac_rr <- c(1)
#'
#' rr_transform <- list(
#'   function(x) {return(c(0.90*(1 - x), (1-0.90)*(1-x), (1-0.90)*x, 0.90*x))}
#' )
#' prob_rr_transform <- c(1)
#'
#' cohorts_max <- 20
#' trial_struc <- "all_plac"
#' safety_prob <- 0
#' sharing_type <- "dynamic"
#' sr_drugs_pos <- 7
#' n_int <- 100
#' n_fin <- 200
#' stage_data <- TRUE
#' cohort_random <- 0.02
#' target_rr <- c(0,0,1)
#' cohort_offset <- 0
#' random_type <- "risk_ratio"
#' sr_first_pos <- FALSE
#'
#' # Vergleich Combo vs Mono
#' Bayes_Sup1 <- matrix(nrow = 1, ncol = 3)
#' Bayes_Sup1[1,] <- c(0.00, 0.90, 1.00)
#' # Vergleich Combo vs Backbone
#' Bayes_Sup2 <- matrix(nrow = 1, ncol = 3)
#' Bayes_Sup2[1,] <- c(0.00, 0.90, 1.00)
#' # Vergleich Mono vs Placebo
#' Bayes_Sup3 <- matrix(nrow = 1, ncol = 3)
#' Bayes_Sup3[1,] <- c(0.00, 0.80, 1.00)
#' Bayes_Sup4 <- matrix(nrow = 1, ncol = 3)
#' Bayes_Sup4[1,] <- c(0.00, 0.80, 1.00)
#' Bayes_Sup <- list(list(Bayes_Sup1, Bayes_Sup2, Bayes_Sup3, Bayes_Sup4),
#'              list(Bayes_Sup1, Bayes_Sup2, Bayes_Sup3, Bayes_Sup4))
#'
#' # Vergleich Combo vs Mono
#' Bayes_Fut1 <- matrix(nrow = 1, ncol = 2)
#' Bayes_Fut1[1,] <- c(0.00, 0.50)
#' # Vergleich Combo vs Backbone
#' Bayes_Fut2 <- matrix(nrow = 1, ncol = 2)
#' Bayes_Fut2[1,] <- c(0.00, 0.50)
#' # Vergleich Mono vs Placebo
#' Bayes_Fut3 <- matrix(nrow = 1, ncol = 2)
#' Bayes_Fut3[1,] <- c(0.00, 0.50)
#' Bayes_Fut4 <- matrix(nrow = 1, ncol = 2)
#' Bayes_Fut4[1,] <- c(0.00, 0.50)
#' Bayes_Fut <- list(list(Bayes_Fut1, Bayes_Fut2, Bayes_Fut3, Bayes_Fut4),
#'                   list(Bayes_Fut1, Bayes_Fut2, Bayes_Fut3, Bayes_Fut4))
#'
#' res_list <- simulate_trial(
#' n_int = n_int, n_fin = n_fin, trial_struc = trial_struc, random_type = random_type,
#' rr_comb = rr_comb, rr_mono = rr_mono, rr_back = rr_back, rr_plac = rr_plac,
#' rr_transform = rr_transform, random = random, prob_comb_rr = prob_comb_rr,
#' prob_mono_rr = prob_mono_rr, prob_back_rr = prob_back_rr, prob_plac_rr = prob_plac_rr,
#' stage_data = stage_data, cohort_random = cohort_random, cohorts_max = cohorts_max,
#' sr_drugs_pos = sr_drugs_pos, target_rr = target_rr, sharing_type = sharing_type,
#' safety_prob = safety_prob, Bayes_Sup = Bayes_Sup, prob_rr_transform = prob_rr_transform,
#' cohort_offset = cohort_offset, Bayes_Fut = Bayes_Fut, sr_first_pos = sr_first_pos
#' )
#'
#' plot_trial(res_list, unit = "n")
#'
#' @export
plot_trial <- function(res_list, unit = "cohort") {

  if (unit == "cohort") {

    # Create Data Frame with several observations per arm
    # Have columns "arm" (Exp1, Contr1, ...), "rr", "Analysis" (how manieth analysis is this?), "N" (up to this point),
    # "Resp" (up to this point), "Decision", further columns for all decision criteria

    "%>%" <- dplyr::"%>%"

    ############## Plot 1 #######

    dat1 <- dplyr::tibble(I = 1:length(res_list$Stage_Data))

    end_n <- function(x, y) {
      ret <- rep(NA, length(x))
      for (i in 1:length(x)) {
        if (is.na(x[i])) {
          ret[i] <- y[i]
        } else {
          ret[i] <- x[i]
        }
      }
      return(ret)
    }

    dat1 <-
      dat1 %>%
      dplyr::mutate(
        Cohort = factor(names(res_list$Stage_Data)),
        Cohort = factor(Cohort, levels = rev(levels(Cohort))),
        # RR_Comb = res_list$Trial_Overview$RR_Comb,
        # RR_Mono = res_list$Trial_Overview$RR_Mono,
        # RR_Back = res_list$Trial_Overview$RR_Back,
        # RR_Plac = res_list$Trial_Overview$RR_Plac,
        Decision_Int = res_list$Trial_Overview$Decision[1,],
        Decision_Fin = res_list$Trial_Overview$Decision[2,],
        Start_N = rep(0, length(Cohort)),
        Interim_N = purrr::map(res_list$Stage_Data, ~ .x$"interim_n_cohort") %>% purrr::flatten_dbl(),
        Final_N = purrr::map(res_list$Stage_Data, ~ .x$"final_n_cohort") %>% purrr::flatten_dbl(),
        End_N = end_n(Final_N, Interim_N),
        Final_Suc = factor(purrr::map(res_list$Stage_Data, ~ .x$"sup_final") %>% purrr::flatten_dbl(), levels = c(0,1)),
        Int_Suc = factor(purrr::map(res_list$Stage_Data, ~ .x$"sup_interim") %>% purrr::flatten_dbl() -
                           purrr::map(res_list$Stage_Data, ~ .x$"fut_interim") %>% purrr::flatten_dbl(),
                         levels = c(-1, 0, 1))
      )
    dat1$Interim_N[dat1$Interim_N == dat1$Final_N] <- NA

    study_design <-
      ggplot2::ggplot(dat1, ggplot2::aes(Decision_Int = Decision_Int, Decision_Fin = Decision_Fin)) +
      ggplot2::geom_segment(ggplot2::aes(x = Cohort, xend = Cohort, y = Start_N, yend = End_N), color = "grey") +
      ggplot2::geom_point(ggplot2::aes(x = Cohort, y = Start_N), size = 2) +
      ggplot2::geom_point(ggplot2::aes(x = Cohort, y = Final_N, fill = Final_Suc), size = 2) +
      ggplot2::scale_fill_manual(values = c("red", "green"), drop = FALSE, guide = FALSE) +
      ggplot2::geom_point(ggplot2::aes(x = Cohort, y = Interim_N, color = Int_Suc), size = 2) +
      ggplot2::scale_color_manual(values = c("red", "orange", "green"), drop = FALSE, guide = FALSE) +
      ggplot2::theme_minimal() +
      ggplot2::xlab("") +
      ggplot2::ylab("N") +
      ggplot2::coord_flip() +
      ggplot2::ggtitle("Overview of Study")

    ply1 <- plotly::ggplotly(study_design) %>% plotly::hide_legend()


    # For future plots:

    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
    }

    CohortColors <-
      stats::setNames(gg_color_hue(length(dat1$Cohort)), levels(dat1$Cohort))

    ########## Plot 2 - Correlation of binary endpoints #########

    bio_comb <- purrr::map(res_list$Stage_Data, ~ .x$Comb$"resp_bio") %>% purrr::flatten_dbl()
    bio_mono <- purrr::map(res_list$Stage_Data, ~ .x$Mono$"resp_bio") %>% purrr::flatten_dbl()
    bio_back <- purrr::map(res_list$Stage_Data, ~ .x$Back$"resp_bio") %>% purrr::flatten_dbl()
    bio_plac <- purrr::map(res_list$Stage_Data, ~ .x$Plac$"resp_bio") %>% purrr::flatten_dbl()
    hist_comb <- purrr::map(res_list$Stage_Data, ~ .x$Comb$"resp_hist") %>% purrr::flatten_dbl()
    hist_mono <- purrr::map(res_list$Stage_Data, ~ .x$Mono$"resp_hist") %>% purrr::flatten_dbl()
    hist_back <- purrr::map(res_list$Stage_Data, ~ .x$Back$"resp_hist") %>% purrr::flatten_dbl()
    hist_plac <- purrr::map(res_list$Stage_Data, ~ .x$Plac$"resp_hist") %>% purrr::flatten_dbl()
    n_comb <- purrr::map(res_list$Stage_Data, ~ .x$Comb$"n") %>% purrr::flatten_dbl()
    n_mono <- purrr::map(res_list$Stage_Data, ~ .x$Mono$"n") %>% purrr::flatten_dbl()
    n_back <- purrr::map(res_list$Stage_Data, ~ .x$Back$"n") %>% purrr::flatten_dbl()
    n_plac <- purrr::map(res_list$Stage_Data, ~ .x$Plac$"n") %>% purrr::flatten_dbl()

    dat2 <-
      dplyr::tibble(
        Bio = c(bio_comb, bio_mono, bio_back, bio_plac),
        Hist = c(hist_comb, hist_mono, hist_back, hist_plac),
        N = c(n_comb, n_mono, n_back, n_plac)
      )

    Bio <- NULL
    Hist <- NULL
    for (i in 1:nrow(dat2)) {
      if (!is.na(dat2$N[i])) {
        B <- numeric(dat2$N[i])
        H <- numeric(dat2$N[i])
        B[0:dat2$Bio[i]] <- 1
        H[0:dat2$Hist[i]] <- 1
        Bio <- c(Bio, B)
        Hist <- c(Hist, H)
      }
    }

    dat2_n <- dplyr::tibble(
      Biomarker_Response = Bio,
      Histology_Response = Hist
    )

    correlation <-
      ggplot2::ggplot(dat2_n, ggplot2::aes(x = Biomarker_Response, y = Histology_Response)) +
      ggplot2::geom_jitter(size = 0.75) +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle("Correlation of Biomarker and Histology Endpoints")

    ply2 <- plotly::ggplotly(correlation)

    ########## Plot 3 - Responders of backbone arms w/r to cohort #########

    # one observation for every arm per length of n-vector
    dat3 <- dplyr::tibble(I = 1:((length(res_list$Stage_Data)*length(res_list$Stage_Data[[1]]$Comb$n))))

    cumsum_help <- function(x) {
      cumsum(ifelse(is.na(x), 0, x))
    }

    dat3 <-
      dat3 %>%
      dplyr::mutate(
        Cohort = forcats::fct_reorder(factor(rep(names(res_list$Stage_Data), each = nrow(dat3)/length(res_list$Stage_Data))), rev(I)),
        N_Cohort = purrr::map(purrr::map(res_list$Stage_Data, ~ .x$Back$"n"), cumsum_help) %>% purrr::flatten_dbl(),
        Resp_Bio = purrr::map(purrr::map(res_list$Stage_Data, ~ .x$Back$"resp_bio"), cumsum_help) %>% purrr::flatten_dbl(),
        Resp_Hist = purrr::map(purrr::map(res_list$Stage_Data, ~ .x$Back$"resp_hist"), cumsum_help) %>% purrr::flatten_dbl(),
        Bio = Resp_Bio / N_Cohort,
        Hist = Resp_Hist / N_Cohort,
        True = rep(purrr::map(res_list$Stage_Data, ~ .x$Back$"rr") %>% purrr::flatten_dbl(), each = nrow(dat3)/length(res_list$Stage_Data))
      ) %>%
      tidyr::gather(key = "Endpoint_Type", value = "RR_Backbone", Bio, Hist, True)


    back_bio_wr_stage <-
      ggplot2::ggplot(dat3, ggplot2::aes(x = N_Cohort, y = RR_Backbone)) +
      ggplot2::geom_line(ggplot2::aes(color = Cohort, linetype = Endpoint_Type)) +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle("Backbone Response Rate") +
      ggplot2::ylim(0, 1) +
      ggplot2::scale_color_manual(values = CohortColors)

    ply3 <- plotly::ggplotly(back_bio_wr_stage)
    ply3$x$layout$annotations[[1]]$text <- ""

    ########## Plot 4 - Responders of placebo arms w/r to cohort #######

    # get only cohorts that have placebo

    have_plac <- sapply(res_list$Stage_Data, function(x) "Plac" %in% names(x))
    if (!any(have_plac)) {
      # empty plot
      df <- data.frame()
      ply4 <- plotly::ggplotly(ggplot2::ggplot(df) + ggplot2::geom_point()) + ggplot2::theme_minimal()
    } else {
      sublist_plac <- res_list$Stage_Data[have_plac]

      # one observation for every arm per length of n-vector
      dat4 <- dplyr::tibble(I = 1:((length(sublist_plac)*length(sublist_plac[[1]]$Plac$n))))

      cumsum_help <- function(x) {
        cumsum(ifelse(is.na(x), 0, x))
      }

      dat4 <-
        dat4 %>%
        dplyr::mutate(
          Cohort = forcats::fct_reorder(factor(rep(names(sublist_plac), each = nrow(dat4)/length(sublist_plac))), rev(I)),
          N_Cohort = purrr::map(purrr::map(sublist_plac, ~ .x$Plac$"n"), cumsum_help) %>% purrr::flatten_dbl(),
          Resp_Bio = purrr::map(purrr::map(sublist_plac, ~ .x$Plac$"resp_bio"), cumsum_help) %>% purrr::flatten_dbl(),
          Resp_Hist = purrr::map(purrr::map(sublist_plac, ~ .x$Plac$"resp_hist"), cumsum_help) %>% purrr::flatten_dbl(),
          Bio = Resp_Bio / N_Cohort,
          Hist = Resp_Hist / N_Cohort,
          True = rep(purrr::map(sublist_plac, ~ .x$Plac$"rr") %>% purrr::flatten_dbl(), each = nrow(dat4)/length(sublist_plac))
        ) %>%
        tidyr::gather(key = "Endpoint_Type", value = "RR_Placebo", Bio, Hist, True)

      if (!all(have_plac)) {
        dat4 <-
          dat4 %>%
          dplyr::mutate(
            Cohort = factor(Cohort, levels = rev(paste0("Cohort", 1:length(have_plac))))
          )
        dat4 <- rbind(dat4,
                      list(I = nrow(dat4) + 1, Cohort = "Cohort1", N_Cohort = 0, Resp_Bio = 0, Resp_Hist = 0, Endpoint_Type = "Bio", RR_Placebo = 0),
                      list(I = nrow(dat4) + 2, Cohort = "Cohort1", N_Cohort = 0, Resp_Bio = 0, Resp_Hist = 0, Endpoint_Type = "Hist", RR_Placebo = 0))
      }

      back_bio_wr_stage <-
        ggplot2::ggplot(dat4, ggplot2::aes(x = N_Cohort, y = RR_Placebo)) +
        ggplot2::geom_line(ggplot2::aes(color = Cohort, linetype = Endpoint_Type)) +
        ggplot2::theme_minimal() +
        ggplot2::ggtitle("Placebo Response Rate") +
        ggplot2::ylim(0, 1) +
        ggplot2::scale_color_manual(values = CohortColors)

      ply4 <- plotly::ggplotly(back_bio_wr_stage)
      ply4$x$layout$annotations[[1]]$text <- ""
    }


    ########## Plot 5 - Responders of combo arms w/r to cohort #########

    # one observation for every arm per length of n-vector
    dat5 <- dplyr::tibble(I = 1:((length(res_list$Stage_Data)*length(res_list$Stage_Data[[1]]$Comb$n))))

    cumsum_help <- function(x) {
      cumsum(ifelse(is.na(x), 0, x))
    }

    dat5 <-
      dat5 %>%
      dplyr::mutate(
        Cohort = forcats::fct_reorder(factor(rep(names(res_list$Stage_Data), each = nrow(dat5)/length(res_list$Stage_Data))), rev(I)),
        N_Cohort = purrr::map(purrr::map(res_list$Stage_Data, ~ .x$Comb$"n"), cumsum_help) %>% purrr::flatten_dbl(),
        Resp_Bio = purrr::map(purrr::map(res_list$Stage_Data, ~ .x$Comb$"resp_bio"), cumsum_help) %>% purrr::flatten_dbl(),
        Resp_Hist = purrr::map(purrr::map(res_list$Stage_Data, ~ .x$Comb$"resp_hist"), cumsum_help) %>% purrr::flatten_dbl(),
        Bio = Resp_Bio / N_Cohort,
        Hist = Resp_Hist / N_Cohort,
        True = rep(purrr::map(res_list$Stage_Data, ~ .x$Comb$"rr") %>% purrr::flatten_dbl(), each = nrow(dat5)/length(res_list$Stage_Data))
      ) %>%
      tidyr::gather(key = "Endpoint_Type", value = "RR_Combo", Bio, Hist, True)


    back_bio_wr_stage <-
      ggplot2::ggplot(dat5, ggplot2::aes(x = N_Cohort, y = RR_Combo)) +
      ggplot2::geom_line(ggplot2::aes(color = Cohort, linetype = Endpoint_Type)) +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle("Combo Response Rate") +
      ggplot2::ylim(0, 1) +
      ggplot2::scale_color_manual(values = CohortColors)

    ply5 <- plotly::ggplotly(back_bio_wr_stage)
    ply5$x$layout$annotations[[1]]$text <- ""


    ########## Plot 6 - Responders of mono arms w/r to cohort #########

    # one observation for every arm per length of n-vector
    dat6 <- dplyr::tibble(I = 1:((length(res_list$Stage_Data)*length(res_list$Stage_Data[[1]]$Comb$n))))

    cumsum_help <- function(x) {
      cumsum(ifelse(is.na(x), 0, x))
    }

    dat6 <-
      dat6 %>%
      dplyr::mutate(
        Cohort = forcats::fct_reorder(factor(rep(names(res_list$Stage_Data), each = nrow(dat6)/length(res_list$Stage_Data))), rev(I)),
        N_Cohort = purrr::map(purrr::map(res_list$Stage_Data, ~ .x$Mono$"n"), cumsum_help) %>% purrr::flatten_dbl(),
        Resp_Bio = purrr::map(purrr::map(res_list$Stage_Data, ~ .x$Mono$"resp_bio"), cumsum_help) %>% purrr::flatten_dbl(),
        Resp_Hist = purrr::map(purrr::map(res_list$Stage_Data, ~ .x$Mono$"resp_hist"), cumsum_help) %>% purrr::flatten_dbl(),
        Bio = Resp_Bio / N_Cohort,
        Hist = Resp_Hist / N_Cohort,
        True = rep(purrr::map(res_list$Stage_Data, ~ .x$Mono$"rr") %>% purrr::flatten_dbl(), each = nrow(dat6)/length(res_list$Stage_Data))
      ) %>%
      tidyr::gather(key = "Endpoint_Type", value = "RR_Mono", Bio, Hist, True)


    back_bio_wr_stage <-
      ggplot2::ggplot(dat6, ggplot2::aes(x = N_Cohort, y = RR_Mono)) +
      ggplot2::geom_line(ggplot2::aes(color = Cohort, linetype = Endpoint_Type)) +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle("Mono Response Rate") +
      ggplot2::ylim(0, 1) +
      ggplot2::scale_color_manual(values = CohortColors)

    ply6 <- plotly::ggplotly(back_bio_wr_stage)
    ply6$x$layout$annotations[[1]]$text <- ""

  } else {

    # Create Data Frame with several observations per arm
    # Have columns "arm" (Exp1, Contr1, ...), "rr", "Analysis" (how manieth analysis is this?), "N" (up to this point),
    # "Resp" (up to this point), "Decision", further columns for all decision criteria

    "%>%" <- dplyr::"%>%"

    ############## Plot 1 #######

    dat1 <- dplyr::tibble(I = 1:length(res_list$Stage_Data))

    end_n <- function(x, y) {
      ret <- rep(NA, length(x))
      for (i in 1:length(x)) {
        if (is.na(x[i])) {
          ret[i] <- y[i]
        } else {
          ret[i] <- x[i]
        }
      }
      return(ret)
    }

    dat1 <-
      dat1 %>%
      dplyr::mutate(
        Cohort = factor(names(res_list$Stage_Data)),
        Cohort = factor(Cohort, levels = rev(levels(Cohort))),
        # RR_Comb = res_list$Trial_Overview$RR_Comb,
        # RR_Mono = res_list$Trial_Overview$RR_Mono,
        # RR_Back = res_list$Trial_Overview$RR_Back,
        # RR_Plac = res_list$Trial_Overview$RR_Plac,
        Decision_Int = res_list$Trial_Overview$Decision[1,],
        Decision_Fin = res_list$Trial_Overview$Decision[2,],
        Start_N = purrr::map(res_list$Stage_Data, ~ .x$"start_n") %>% purrr::flatten_dbl(),
        Interim_N = purrr::map(res_list$Stage_Data, ~ .x$"interim_n") %>% purrr::flatten_dbl(),
        Final_N = purrr::map(res_list$Stage_Data, ~ .x$"final_n") %>% purrr::flatten_dbl(),
        End_N = end_n(Final_N, Interim_N),
        Final_Suc = factor(purrr::map(res_list$Stage_Data, ~ .x$"sup_final") %>% purrr::flatten_dbl(), levels = c(0,1)),
        Int_Suc = factor(purrr::map(res_list$Stage_Data, ~ .x$"sup_interim") %>% purrr::flatten_dbl() -
                           purrr::map(res_list$Stage_Data, ~ .x$"fut_interim") %>% purrr::flatten_dbl(),
                         levels = c(-1, 0, 1))
      )
    dat1$Interim_N[dat1$Interim_N == dat1$Final_N] <- NA

    study_design <-
      ggplot2::ggplot(dat1, ggplot2::aes(Decision_Int = Decision_Int, Decision_Fin = Decision_Fin)) +
      ggplot2::geom_segment(ggplot2::aes(x = Cohort, xend = Cohort, y = Start_N, yend = End_N), color = "grey") +
      ggplot2::geom_point(ggplot2::aes(x = Cohort, y = Start_N), size = 2) +
      ggplot2::geom_point(ggplot2::aes(x = Cohort, y = Final_N, fill = Final_Suc), size = 2) +
      ggplot2::scale_fill_manual(values = c("red", "green"), drop = FALSE, guide = FALSE) +
      ggplot2::geom_point(ggplot2::aes(x = Cohort, y = Interim_N, color = Int_Suc), size = 2) +
      ggplot2::scale_color_manual(values = c("red", "orange", "lightgreen"), drop = FALSE, guide = FALSE) +
      ggplot2::theme_minimal() +
      ggplot2::xlab("") +
      ggplot2::ylab("N") +
      ggplot2::coord_flip() +
      ggplot2::ggtitle("Overview of Study")

    ply1 <- plotly::ggplotly(study_design)

    # For future plots:

    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
    }

    CohortColors <-
      stats::setNames(gg_color_hue(length(dat1$Cohort)), levels(dat1$Cohort))

    ########## Plot 2 - Correlation of binary endpoints #########

    bio_comb <- purrr::map(res_list$Stage_Data, ~ .x$Comb$"resp_bio") %>% purrr::flatten_dbl()
    bio_mono <- purrr::map(res_list$Stage_Data, ~ .x$Mono$"resp_bio") %>% purrr::flatten_dbl()
    bio_back <- purrr::map(res_list$Stage_Data, ~ .x$Back$"resp_bio") %>% purrr::flatten_dbl()
    bio_plac <- purrr::map(res_list$Stage_Data, ~ .x$Plac$"resp_bio") %>% purrr::flatten_dbl()
    hist_comb <- purrr::map(res_list$Stage_Data, ~ .x$Comb$"resp_hist") %>% purrr::flatten_dbl()
    hist_mono <- purrr::map(res_list$Stage_Data, ~ .x$Mono$"resp_hist") %>% purrr::flatten_dbl()
    hist_back <- purrr::map(res_list$Stage_Data, ~ .x$Back$"resp_hist") %>% purrr::flatten_dbl()
    hist_plac <- purrr::map(res_list$Stage_Data, ~ .x$Plac$"resp_hist") %>% purrr::flatten_dbl()
    n_comb <- purrr::map(res_list$Stage_Data, ~ .x$Comb$"n") %>% purrr::flatten_dbl()
    n_mono <- purrr::map(res_list$Stage_Data, ~ .x$Mono$"n") %>% purrr::flatten_dbl()
    n_back <- purrr::map(res_list$Stage_Data, ~ .x$Back$"n") %>% purrr::flatten_dbl()
    n_plac <- purrr::map(res_list$Stage_Data, ~ .x$Plac$"n") %>% purrr::flatten_dbl()

    dat2 <-
      dplyr::tibble(
        Bio = c(bio_comb, bio_mono, bio_back, bio_plac),
        Hist = c(hist_comb, hist_mono, hist_back, hist_plac),
        N = c(n_comb, n_mono, n_back, n_plac)
      )

    Bio <- NULL
    Hist <- NULL
    for (i in 1:nrow(dat2)) {
      if (!is.na(dat2$N[i])) {
        B <- numeric(dat2$N[i])
        H <- numeric(dat2$N[i])
        B[0:dat2$Bio[i]] <- 1
        H[0:dat2$Hist[i]] <- 1
        Bio <- c(Bio, B)
        Hist <- c(Hist, H)
      }
    }

    dat2_n <- dplyr::tibble(
      Biomarker_Response = Bio,
      Histology_Response = Hist
    )

    correlation <-
      ggplot2::ggplot(dat2_n, ggplot2::aes(x = Biomarker_Response, y = Histology_Response)) +
      ggplot2::geom_jitter(size = 0.75) +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle("Correlation of Biomarker and Histology Endpoints")

    ply2 <- plotly::ggplotly(correlation)

    ########## Plot 3 - Responders of backbone arms w/r to cohort #########

    # one observation for every arm per length of n-vector
    dat3 <- dplyr::tibble(I = 1:((length(res_list$Stage_Data)*length(res_list$Stage_Data[[1]]$Comb$n))))

    cumsum_help <- function(x) {
      cumsum(ifelse(is.na(x), 0, x))
    }

    dat3 <-
      dat3 %>%
      dplyr::mutate(
        Cohort = forcats::fct_reorder(factor(rep(names(res_list$Stage_Data), each = nrow(dat3)/length(res_list$Stage_Data))), rev(I)),
        N_Total = rep(res_list$Trial_Overview$Total_N_Vector, length(res_list$Stage_Data)),
        N_Cohort = purrr::map(purrr::map(res_list$Stage_Data, ~ .x$Back$"n"), cumsum_help) %>% purrr::flatten_dbl(),
        Resp_Bio = purrr::map(purrr::map(res_list$Stage_Data, ~ .x$Back$"resp_bio"), cumsum_help) %>% purrr::flatten_dbl(),
        Resp_Hist = purrr::map(purrr::map(res_list$Stage_Data, ~ .x$Back$"resp_hist"), cumsum_help) %>% purrr::flatten_dbl(),
        Bio = Resp_Bio / N_Cohort,
        Hist = Resp_Hist / N_Cohort,
        True = rep(purrr::map(res_list$Stage_Data, ~ .x$Back$"rr") %>% purrr::flatten_dbl(), each = nrow(dat3)/length(res_list$Stage_Data))
      ) %>%
      tidyr::gather(key = "Endpoint_Type", value = "RR_Backbone", Bio, Hist, True)


    back_bio_wr_stage <-
      ggplot2::ggplot(dat3, ggplot2::aes(x = N_Total, y = RR_Backbone)) +
      ggplot2::geom_line(ggplot2::aes(color = Cohort, linetype = Endpoint_Type)) +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle("Backbone Response Rate") +
      ggplot2::ylim(0, 1) +
      ggplot2::scale_color_manual(values = CohortColors)

    ply3 <- plotly::ggplotly(back_bio_wr_stage)
    ply3$x$layout$annotations[[1]]$text <- ""

    ########## Plot 4 - Responders of placebo arms w/r to cohort #######

    # get only cohorts that have placebo

    have_plac <- sapply(res_list$Stage_Data, function(x) "Plac" %in% names(x))
    if (!any(have_plac)) {
      # empty plot
      df <- data.frame()
      ply4 <- plotly::ggplotly(ggplot2::ggplot(df) + ggplot2::geom_point()) + ggplot2::theme_minimal()
    } else {
      sublist_plac <- res_list$Stage_Data[have_plac]

      # one observation for every arm per length of n-vector
      dat4 <- dplyr::tibble(I = 1:((length(sublist_plac)*length(sublist_plac[[1]]$Plac$n))))

      cumsum_help <- function(x) {
        cumsum(ifelse(is.na(x), 0, x))
      }

      dat4 <-
        dat4 %>%
        dplyr::mutate(
          Cohort = forcats::fct_reorder(factor(rep(names(sublist_plac), each = nrow(dat4)/length(sublist_plac))), rev(I)),
          N_Cohort = purrr::map(purrr::map(sublist_plac, ~ .x$Plac$"n"), cumsum_help) %>% purrr::flatten_dbl(),
          N_Total = rep(res_list$Trial_Overview$Total_N_Vector, length(sublist_plac)),
          Resp_Bio = purrr::map(purrr::map(sublist_plac, ~ .x$Plac$"resp_bio"), cumsum_help) %>% purrr::flatten_dbl(),
          Resp_Hist = purrr::map(purrr::map(sublist_plac, ~ .x$Plac$"resp_hist"), cumsum_help) %>% purrr::flatten_dbl(),
          Bio = Resp_Bio / N_Cohort,
          Hist = Resp_Hist / N_Cohort,
          True = rep(purrr::map(sublist_plac, ~ .x$Plac$"rr") %>% purrr::flatten_dbl(), each = nrow(dat4)/length(sublist_plac))
        ) %>%
        tidyr::gather(key = "Endpoint_Type", value = "RR_Placebo", Bio, Hist, True)

      if (!all(have_plac)) {
        dat4 <-
          dat4 %>%
          dplyr::mutate(
            Cohort = factor(Cohort, levels = rev(paste0("Cohort", 1:length(have_plac))))
          )
        dat4 <- rbind(dat4,
                      list(I = nrow(dat4) + 1, Cohort = "Cohort1", N_Cohort = 0, N_Total = 0, Resp_Bio = 0, Resp_Hist = 0, Endpoint_Type = "Bio", RR_Placebo = 0),
                      list(I = nrow(dat4) + 2, Cohort = "Cohort1", N_Cohort = 0, N_Total = 0, Resp_Bio = 0, Resp_Hist = 0, Endpoint_Type = "Hist", RR_Placebo = 0))
      }

      back_bio_wr_stage <-
        ggplot2::ggplot(dat4, ggplot2::aes(x = N_Total, y = RR_Placebo)) +
        ggplot2::geom_line(ggplot2::aes(color = Cohort, linetype = Endpoint_Type)) +
        ggplot2::theme_minimal() +
        ggplot2::ggtitle("Placebo Response Rate") +
        ggplot2::ylim(0, 1) +
        ggplot2::scale_color_manual(values = CohortColors)

      ply4 <- plotly::ggplotly(back_bio_wr_stage)
      ply4$x$layout$annotations[[1]]$text <- ""
    }


    ########## Plot 5 - Responders of combo arms w/r to cohort #########

    # one observation for every arm per length of n-vector
    dat5 <- dplyr::tibble(I = 1:((length(res_list$Stage_Data)*length(res_list$Stage_Data[[1]]$Comb$n))))

    cumsum_help <- function(x) {
      cumsum(ifelse(is.na(x), 0, x))
    }

    dat5 <-
      dat5 %>%
      dplyr::mutate(
        Cohort = forcats::fct_reorder(factor(rep(names(res_list$Stage_Data), each = nrow(dat5)/length(res_list$Stage_Data))), rev(I)),
        N_Cohort = purrr::map(purrr::map(res_list$Stage_Data, ~ .x$Comb$"n"), cumsum_help) %>% purrr::flatten_dbl(),
        N_Total = rep(res_list$Trial_Overview$Total_N_Vector, length(res_list$Stage_Data)),
        Resp_Bio = purrr::map(purrr::map(res_list$Stage_Data, ~ .x$Comb$"resp_bio"), cumsum_help) %>% purrr::flatten_dbl(),
        Resp_Hist = purrr::map(purrr::map(res_list$Stage_Data, ~ .x$Comb$"resp_hist"), cumsum_help) %>% purrr::flatten_dbl(),
        Bio = Resp_Bio / N_Cohort,
        Hist = Resp_Hist / N_Cohort,
        True = rep(purrr::map(res_list$Stage_Data, ~ .x$Comb$"rr") %>% purrr::flatten_dbl(), each = nrow(dat5)/length(res_list$Stage_Data))
      ) %>%
      tidyr::gather(key = "Endpoint_Type", value = "RR_Combo", Bio, Hist, True)


    back_bio_wr_stage <-
      ggplot2::ggplot(dat5, ggplot2::aes(x = N_Total, y = RR_Combo)) +
      ggplot2::geom_line(ggplot2::aes(color = Cohort, linetype = Endpoint_Type)) +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle("Combo Response Rate") +
      ggplot2::ylim(0, 1) +
      ggplot2::scale_color_manual(values = CohortColors)

    ply5 <- plotly::ggplotly(back_bio_wr_stage)
    ply5$x$layout$annotations[[1]]$text <- ""


    ########## Plot 6 - Responders of mono arms w/r to cohort #########

    # one observation for every arm per length of n-vector
    dat6 <- dplyr::tibble(I = 1:((length(res_list$Stage_Data)*length(res_list$Stage_Data[[1]]$Comb$n))))

    cumsum_help <- function(x) {
      cumsum(ifelse(is.na(x), 0, x))
    }

    dat6 <-
      dat6 %>%
      dplyr::mutate(
        Cohort = forcats::fct_reorder(factor(rep(names(res_list$Stage_Data), each = nrow(dat6)/length(res_list$Stage_Data))), rev(I)),
        N_Cohort = purrr::map(purrr::map(res_list$Stage_Data, ~ .x$Mono$"n"), cumsum_help) %>% purrr::flatten_dbl(),
        N_Total = rep(res_list$Trial_Overview$Total_N_Vector, length(res_list$Stage_Data)),
        Resp_Bio = purrr::map(purrr::map(res_list$Stage_Data, ~ .x$Mono$"resp_bio"), cumsum_help) %>% purrr::flatten_dbl(),
        Resp_Hist = purrr::map(purrr::map(res_list$Stage_Data, ~ .x$Mono$"resp_hist"), cumsum_help) %>% purrr::flatten_dbl(),
        Bio = Resp_Bio / N_Cohort,
        Hist = Resp_Hist / N_Cohort,
        True = rep(purrr::map(res_list$Stage_Data, ~ .x$Mono$"rr") %>% purrr::flatten_dbl(), each = nrow(dat6)/length(res_list$Stage_Data))
      ) %>%
      tidyr::gather(key = "Endpoint_Type", value = "RR_Mono", Bio, Hist, True)


    back_bio_wr_stage <-
      ggplot2::ggplot(dat6, ggplot2::aes(x = N_Total, y = RR_Mono)) +
      ggplot2::geom_line(ggplot2::aes(color = Cohort, linetype = Endpoint_Type)) +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle("Mono Response Rate") +
      ggplot2::ylim(0, 1) +
      ggplot2::scale_color_manual(values = CohortColors)

    ply6 <- plotly::ggplotly(back_bio_wr_stage)
    ply6$x$layout$annotations[[1]]$text <- ""

  }

 ##### Wrap up #####

  p <- plotly::subplot(ply1, ply2, ply5, ply6, ply3, ply4, nrows = 3, titleX = TRUE, titleY = TRUE, margin = 0.05)
  #p$x$layout$title$text <- "a) Study Overview, b) Correlation Bio/Hist,
  #c) CombRR/Cohort, d) MonoRR/Cohort, e) BackRR/Cohort, f) PlacRR/Cohort"
  p$x$layout$title$text <- ""

  nam <- purrr::map(p$x$data, ~ .x$"name") %>% purrr::flatten_chr()
  ind_leg <- NULL
  ind_names <- NULL
  for (i in 1:length(nam)) {
    if (substr(nam[i], 1, 2) == "(C") {
      if (!nam[i] %in% ind_names) {
        ind_leg <- c(ind_leg, i)
        ind_names <- c(ind_names, nam[i])
      }
    }
  }

  for(i in setdiff(1:length(nam), ind_leg)) {
    p$x$data[[i]]$showlegend <- FALSE
  }

  p


}
