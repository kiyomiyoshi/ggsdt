#' @export
fit_ggsdt <- function(nR_S1, nR_S2, add_constant = TRUE) {

    if (add_constant) {
        nR_S1 <- nR_S1 + (1 / length(nR_S1))
        nR_S2 <- nR_S2 + (1 / length(nR_S2))
    }

    n_ratings <- length(nR_S1) / 2
    ratingFAR <- 1 - cumsum(nR_S1) / sum(nR_S1)
    ratingHR <-  1 - cumsum(nR_S2) / sum(nR_S2)

    # set up initial guess at parameter values
    alp2 <- 1
    bet <- 2
    d1 <- gnorm::qgnorm(ratingHR[n_ratings], alpha = alp2, beta = bet) -
        gnorm::qgnorm(ratingFAR[n_ratings], alpha = alp2, beta = bet)
    c1 <- -1 * gnorm::qgnorm(ratingFAR, alpha = alp2, beta = bet)
    cri <- c1[1:(2 * n_ratings - 1)]
    guess <- c(d1, alp2, bet, cri) # can be modified for better model convergence

    # model fitting
    params <- list("n_ratings" = n_ratings, "nR_S1" = nR_S1, "nR_S2" = nR_S2)

    fit <- suppressWarnings(stats::optim(par = guess, fit_ggsdt_logL, gr = NULL, method = "BFGS", parameters = params,
                                  lower = c(0, 0, 0, rep(-Inf, 2 * n_ratings - 1)),
                                  control = list("maxit" = 100000,
                                                 "parscale" = c(1, 0.3, 0.3, rep(0.1, 2 * n_ratings - 1)))))

    m2 <-    fit$par[1]
    alp2 <- fit$par[2]
    bet <-  fit$par[3]
    logL <- -fit$value
    sd1 <-  sqrt((1^2 * gamma(3 / bet)) / gamma(1 / bet))
    sd2 <-  sqrt((alp2^2 * gamma(3 / bet)) / gamma(1 / bet))
    kurt <- (gamma(5 / bet) * gamma(1 / bet)) / gamma(3 / bet)^2 - 3

    est <- data.frame(mu2 = m2, alpha2 = alp2, beta = bet, LogLike = logL,
                      sigma1 = sd1, sigma2 = sd2, kurtosis = kurt)
    cri <- data.frame(matrix(vector(), 0, 2 * n_ratings - 1))

    for (i in 1:(2 * n_ratings - 1)) {
        cri[1, i] <- fit$par[3 + i]
    }

    est <- cbind(est, cri)

    return(est)

}


fit_ggsdt_logL <- function(x, parameters) {

    d1 <- x[1]
    alp2 <- x[2]
    bet <- x[3]
    cri <- x[4:(2 * parameters$n_ratings + 2)]

    nR_S1 <- parameters$nR_S1
    nR_S2 <- parameters$nR_S2
    n_ratings <- parameters$n_ratings

    far <- gnorm::pgnorm(cri, alpha = 1, beta = bet, mu = 0)
    hr <- gnorm::pgnorm(cri - d1, alpha = alp2, beta = bet, mu = 0)

    s1_exp <- sum(nR_S1) * diff(far)
    s2_exp <- sum(nR_S2) * diff(hr)

    s1_exp_number <- c(sum(nR_S1) * far[1], s1_exp, sum(nR_S1) - sum(nR_S1) * far[1] - sum(s1_exp))
    s2_exp_number <- c(sum(nR_S2) * hr[1], s2_exp, sum(nR_S2) - sum(nR_S2) * hr[1] - sum(s2_exp))

    logL <- sum(nR_S2 * log(s2_exp_number/sum(nR_S2)) + nR_S1 * log(s1_exp_number/sum(nR_S1)))

    if (is.nan(logL)) {
        logL <- -Inf
    }

    logL <- -logL

    return(logL)

}


#' @export
ggdistr <- function(mu2, alpha2, beta) {

    ggplot2::ggplot() + ggplot2::theme_classic() +
        ggplot2::geom_function(fun = gnorm::dgnorm, size = 0.8,
                              args = list(mu = 0, alpha = 1, beta = beta)) +
        ggplot2::geom_function(fun = gnorm::dgnorm, size = 0.8,
                              args = list(mu = mu2, alpha = alpha2, beta = beta)) +
        ggplot2::annotate("text", x =  2.6, y = 0.3, parse = F, label = "S2", size = 5) +
        ggplot2::annotate("text", x = -1.5, y = 0.3, parse = F, label = "S1", size = 5) +
        ggplot2::scale_x_continuous(breaks = seq(0, mu2, mu2),
                           labels = c(0, round(mu2, digits = 3)), limits = c(-5, 5), expand = c(0, 0)) +
        ggplot2::scale_y_continuous(breaks = NULL,
                           limits = c(0, 0.75), expand = c(0, 0)) +
        ggplot2::xlab("") + ggeasy::easy_remove_y_axis()

}


#' @export
ggroc1 <- function(mu2, alpha2, beta) {

    dat <- c()

    for (i in seq(-10, 10, by = 0.02)) {
        dat <- rbind(dat, c(1 - gnorm::pgnorm(i, mu = 0, alpha = 1, beta = beta),
                            1 - gnorm::pgnorm(i, mu = mu2, alpha = alpha2, beta = beta)))
    }

    dat <- as.data.frame(dat)
    colnames(dat) <- c("far", "hr")

    ggplot2::ggplot(dat) + ggplot2::theme_classic() +
        ggplot2::geom_line(ggplot2::aes(x = far, y = hr)) +
        ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
        ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
        ggplot2::xlab("Hit rate") + ggplot2::ylab("False alarm rate") +
        ggplot2::coord_fixed(ratio = 1)

}


#' @export
ggzroc1 <- function(mu2, alpha2, beta) {

    dat <- c()

    for (i in seq(-10, 10, by = 0.02)) {
        dat <- rbind(dat, c(1 - gnorm::pgnorm(i, mu = 0, alpha = 1, beta = beta),
                            1 - gnorm::pgnorm(i, mu = mu2, alpha = alpha2, beta = beta)))
    }

    dat <- as.data.frame(dat)

    ggplot2::ggplot(dat) + ggplot2::theme_classic() +
        ggplot2::geom_line(ggplot2::aes(x = stats::qnorm(far), y = stats::qnorm(hr))) +
        ggplot2::scale_x_continuous(limits = c(-2.5, 2.5), breaks = seq(-2, 2, by = 1)) +
        ggplot2::scale_y_continuous(limits = c(-2.5, 2.5), breaks = seq(-2, 2, by = 1)) +
        ggplot2::xlab("z(hit rate)") + ggplot2::ylab("z(false alarm rate)") +
        ggplot2::coord_fixed(ratio = 1)

}
