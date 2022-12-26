#' @export
fit_ggsdt <- function(nR_S1, nR_S2, add_constant = TRUE) {

    if (add_constant) {
        nR_S1 <- nR_S1 + (1 / length(nR_S1))
        nR_S2 <- nR_S2 + (1 / length(nR_S2))
    }

    n_ratings <- length(nR_S1) / 2
    far <- 1 - cumsum(nR_S1) / sum(nR_S1)
    hr <-  1 - cumsum(nR_S2) / sum(nR_S2)

    # set up initial guess for parameter values
    alp2 <- 1
    bet <-  2
    m2 <-   gnorm::qgnorm(hr[n_ratings], alpha = alp2, beta = bet) -
            gnorm::qgnorm(far[n_ratings], alpha = alp2, beta = bet)
    cri <-  -1 * gnorm::qgnorm(far, alpha = alp2, beta = bet)
    cri <-  cri[1:(2 * n_ratings - 1)]

    guess <- c(m2, alp2, bet, cri) # can be modified for model convergence

    # model fitting
    params <- list("n_ratings" = n_ratings, "nR_S1" = nR_S1, "nR_S2" = nR_S2)

    fit <- suppressWarnings(stats::optim(par = guess, fn = fit_ggsdt_ll, gr = NULL, method = "L-BFGS-B", parameters = params,
                                  lower = c(0, 0, 0, rep(-Inf, 2 * n_ratings - 1)),
                                  control = list("maxit" = 100000,
                                                 "parscale" = c(1, 0.3, 0.3, rep(0.1, 2 * n_ratings - 1)))))

    m2 <-   fit$par[1]
    alp2 <- fit$par[2]
    bet <-  fit$par[3]
    ll <-   -fit$value
    sd1 <-  sqrt((1^2 * gamma(3 / bet)) / gamma(1 / bet))
    sd2 <-  sqrt((alp2^2 * gamma(3 / bet)) / gamma(1 / bet))
    kurt <- (gamma(5 / bet) * gamma(1 / bet)) / gamma(3 / bet)^2 - 3

    est <- data.frame(mu2 = m2, alpha2 = alp2, beta = bet, loglike = ll,
                      sigma1 = sd1, sigma2 = sd2, kurtosis = kurt)

    for (i in 1:(2 * n_ratings - 1)) {
        new <- fit$par[3 + i]
        est[, ncol(est) + 1] <- new
        colnames(est)[ncol(est)] <- paste0("c", i)
    }

    return(est)

}


fit_ggsdt_ll <- function(x, parameters) {

    m2 <-   x[1]
    alp2 <- x[2]
    bet <-  x[3]
    cri <-  x[4:(2 * parameters$n_ratings + 2)]

    nR_S1 <- parameters$nR_S1
    nR_S2 <- parameters$nR_S2

    exp_far <- gnorm::pgnorm(cri, alpha = 1, beta = bet, mu = 0)
    exp_hr <- gnorm::pgnorm(cri - m2, alpha = alp2, beta = bet, mu = 0)

    exp_s1 <- sum(nR_S1) * diff(exp_far)
    exp_s2 <- sum(nR_S2) * diff(exp_hr)

    exp_fas <-  c(sum(nR_S1) * exp_far[1], exp_s1, sum(nR_S1) - sum(nR_S1) * exp_far[1] - sum(exp_s1))
    exp_hits <- c(sum(nR_S2) * exp_hr[1],  exp_s2, sum(nR_S2) - sum(nR_S2) * exp_hr[1]  - sum(exp_s2))

    ll <- sum(nR_S2 * log(exp_hits / sum(nR_S2)) + nR_S1 * log(exp_fas / sum(nR_S1)))

    if (is.nan(ll)) {
        ll <- -Inf
    }

    ll <- -ll

    return(ll)

}


#' @export
ggdistr <- function(mu2, alpha2, beta) {

    ggplot2::ggplot() + ggplot2::theme_classic() +
        ggplot2::geom_function(fun = gnorm::dgnorm, size = 0.8,
                              args = list(mu = 0, alpha = 1, beta = beta)) +
        ggplot2::geom_function(fun = gnorm::dgnorm, size = 0.8,
                              args = list(mu = mu2, alpha = alpha2, beta = beta)) +
        ggplot2::annotate("text", x =  2.6, y = 0.3, parse = FALSE, label = "S2", size = 5) +
        ggplot2::annotate("text", x = -1.5, y = 0.3, parse = FALSE, label = "S1", size = 5) +
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
    colnames(dat) <- c("far", "hr")

    ggplot2::ggplot(dat) + ggplot2::theme_classic() +
        ggplot2::geom_line(ggplot2::aes(x = stats::qnorm(far), y = stats::qnorm(hr))) +
        ggplot2::scale_x_continuous(limits = c(-2.5, 2.5), breaks = seq(-2, 2, by = 1)) +
        ggplot2::scale_y_continuous(limits = c(-2.5, 2.5), breaks = seq(-2, 2, by = 1)) +
        ggplot2::xlab("z(hit rate)") + ggplot2::ylab("z(false alarm rate)") +
        ggplot2::coord_fixed(ratio = 1)

}
