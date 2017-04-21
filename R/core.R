#' Approximate scores for ranks.
#'
#' This function modifies the \code{multrnks} function in the \code{sensitivitymw} package by also providing the exact scores. The scores are also normalized so that the maximum is 1.
#'
#' @param rk a vector of ranks
#' @inheritParams sen
#'
#' @author Paul Rosenbaum, Qingyuan Zhao
#' @details
#' Exact and approximate rank scores yield similar bounds on P-values when sample size is large.
#' The exact rank scores involve very large combinatorial coefficiences when the same size is very large, whereas the nearly equivalent approximate scores do not.
#'
#' @export
#'
multrnks <- function (rk, mm, score.method = c("approximate", "exact"))
{
    score.method <- match.arg(score.method, c("approximate", "exact"))
    n <- length(rk)
    pk <- rk/n
    q <- rep(0, n)
    m1 <- mm[2]
    m2 <- mm[3]
    m <- mm[1]
    if (score.method == "approximate") {
        for (l in m1:m2) {
            q <- q + (l * choose(m, l) * (pk^(l - 1)) * ((1 - pk)^(m - l)))
        }
    } else if (score.method == "exact") {
        for (l in m1:m2){
            q <- q + choose(rk-1,l-1) * choose(n-rk,m-l)
        }
    } else {
        stop("score.method must either be approximate or exact.")
    }

    q / max(q)
}

#' Sensivity analysis with signed score test
#'
#' This function implements Rosenbaum's sensitivity analysis for pair-matched observational study with general signed score test. It is faster and more flexible than the \code{psens} function in the package \code{rbounds}.
#'
#' @param d a vector of treatment-minus-control differences
#' @param mm a vector (m, munder, mover) or a matrix, each column a vector (m, munder, mover) that indicates the U-statistic.s NULL means Wilcoxon's signed rank test.
#' @param gamma a vector of sensitivity parameters (must be >= 1).
#' @param alternative report p-value corresponds to the maximum ("upper") or minimum ("lower") bound
#' @param approx.method how to compute the $p$-value upper bound? either "normal" approximation or random "permutations".
#' @param score.method either approximate score or exact score
#' @param tau a scalar, null hypothesis is the additive effect is \code{tau} (default 0)
#' @param num.perms number of Monte-Carlo simulations used to compute the sensivitiy value, if \code{approx.method} is "permutations".
#'
#' @importFrom stats pnorm
#'
#' @return A list
#' \describe{
#' \item{p.value}{p-values corresponding to each entry of \code{gamma}}
#' \item{p.value2}{two sided p-values}
#' \item{gamma.hat}{estimate of design sensitivity}
#' \item{T}{test statistic}
#' \item{E}{Means of the test statistic under sensivity \code{gamma}}
#' \item{V}{Variances of the test statistic under sensitivity \code{gamma}}
#' \item{eff.size}{Effect size of T compared to E and V}
#' \item{E.gamma1}{Expectation of T under null at Gamma = 1}
#' }
#'
#' @references
#' \itemize{
#' \item{Rosenbaum, Paul R. \emph{Observational Studies}. Springer New York, 2002.}
#' \item{Rosenbaum, P. R. (2011). A New u-Statistic with Superior Design Sensitivity in Matched Observational Studies. \emph{Biometrics}, 67(3), 1017-1027.}
#' }
#'
#' @author Paul Rosenbaum, Qingyuan Zhao
#' @export
#'
#' @examples
#'
#' require(CrossScreening)
#' data(lead)
#' d.lead <- lead$exposed[-21] - lead$control[-21]
#' sen(d.lead, gamma = c(1, 2, 3, 4, 5, 6))
#'
sen <- function(d, mm = NULL, gamma = 1, alternative = c("greater", "less"), approx.method = c("normal", "permutation"), score.method = c("approximate", "exact"), tau = 0, num.perms = 10000) {

    score.method <- match.arg(score.method, c("approximate", "exact"))
    approx.method <- match.arg(approx.method, c("normal", "permutation"))
    alternative <- match.arg(alternative, c("greater", "less"))
    d <- d - tau
    if (alternative == "less") {
        d <- - d
    }

    if (is.null(mm)) {
        ## compute score based on the rank of abs(d)
        qs <- rank(abs(d))
    } else {
        mm <- as.matrix(mm)
        stopifnot(nrow(mm) == 3)
        if (ncol(mm) == 1) {
            qs <- multrnks(rank(abs(d)), mm, score.method)
        } else {
            sen.result <- apply(mm, 2, function(mmm) sen(d, mmm, gamma, score.method = score.method))
            Map.combine <- function(...) {
                Map(cbind, ...)
            }
            sen.result <- do.call(Map.combine, sen.result)
            mm.names <- paste0("(", apply(mm, 2, paste, collapse = ","), ")")
            for (i in 1:length(sen.result)) {
                colnames(sen.result[[i]]) <- mm.names
            }
            return(sen.result)
        }
    }

    qs <- qs * (abs(d) > 0)
    sg <- as.numeric(d > 0) # sign
    Ts <- sum(qs * sg) # test statistic
    ETGamma1 <- sum(qs)/2 # expectation of the statistic when Gamma = 1

    if (approx.method == "permutation") {
        p.value <- rep(0, length(gamma))
        p.value2 <- rep(0, length(gamma))
        for (i in 1:length(gamma)) {
            s <- replicate(num.perms, rbinom(length(qs), 1, prob = gamma2kappa(gamma[i])))
            p.value[i] <- mean(qs %*% s >= Ts)
            p.value2[i] <- pmin(2 * p.value[i], 2 * mean(qs %*% s >= Ts), 1)
        }
        ET <- NULL
        VT <- NULL
        dev <- NULL
        gammahat <- NULL
        names(p.value) <- gamma
        names(p.value2) <- gamma
    } else {
        kappa <- gamma / (1 + gamma)
        VT <- sum(qs*qs)*kappa*(1-kappa) # variance of the bounding statistic
        ET <- sum(qs) * kappa
        dev <- (Ts - ET) / sqrt(VT)
        p.value <- 1 - pnorm(dev)
        kappahat <- Ts / sum(qs)
        gammahat <- max(1, kappahat / (1 - kappahat))

        p.value2 <- 2 * pmin(p.value, 1 - pnorm((sum(qs) - Ts - ET)/sqrt(VT)), 1/2)
        ## devc <- pmax(0, dev)

        names(p.value) <- gamma
        names(p.value2) <- gamma
        names(ET) <- gamma
        names(VT) <- gamma
        names(dev) <- gamma
    }

    return(list(p.value=p.value,
                p.value2=p.value2,
                gamma.hat=gammahat,
                T=Ts,
                E=ET,
                V=VT,
                eff.size=dev,
                E.gamma1=ETGamma1))
}

#' Point estimate and confidence interval for sensitivity analysis
#'
#' @inheritParams sen
#' @param mm a vector (m, munder, mover) that indicates the U-statistic. Does not support matrix \code{mm} in this function.
#' @param alpha significance level for the outer confidence interval
#' @param alpha.up upper-tail probability of the confidence interval
#' @param alpha.low lower-tail probability of the confidence interval
#'
#' @return a list
#' \describe{
#' \item{point.estimate}{An interval of point estimates allowing for a bias of gamma in treatment assignment.}
#' \item{ci}{An confidence interval allowing for a bias of gamma in treatment assignment.}
#' }
#'
#' @details See the \code{senmwCI} function in the \code{sensitivitymw} package.
#'
#' @author Qingyuan Zhao
#'
#' @export
#'
#' @examples
#' data(lead)
#' d.lead <- lead$exposed[-21] - lead$control[-21]
#' sen.ci(d.lead, gamma = c(1, 2), alpha.up = 0, alpha.low = 0.05)
#'
sen.ci <- function(d, mm = c(2, 2, 2), gamma = 1, alpha = 0.05, alpha.up = alpha/2, alpha.low = alpha/2, score.method = c("approximate", "exact")) {

    mm <- as.matrix(mm)
    stopifnot(nrow(mm) == 3 & ncol(mm) == 1)

    score.method <- match.arg(score.method, c("approximate", "exact"))

    inner.ci.fun <- function(tau, alternative, j) {
        sen(d, mm, gamma, alternative = alternative, score.method = score.method, tau = tau)$eff.size[j]
    }
    inner.ci.low <- sapply(1:length(gamma), function(j) uniroot(inner.ci.fun, range(d), alternative = "greater", j = j)$root)
    inner.ci.up <- sapply(1:length(gamma), function(j) uniroot(inner.ci.fun, range(d), alternative = "less", j = j)$root)

    outer.ci.fun <- function(tau, alternative, j, alpha) {
        sen(d, mm, gamma, alternative = alternative, score.method = score.method, tau = tau)$p.value[j] - alpha
    }
    outer.ci.low <- sapply(1:length(gamma), function(j) tryCatch(uniroot(outer.ci.fun, range(d) + 100 * max(gamma) * sd(d) * c(-1, 1), alternative = "greater", j = j, alpha = alpha.low)$root, error = function(e) -Inf))
    outer.ci.up <- sapply(1:length(gamma), function(j) tryCatch(uniroot(outer.ci.fun, range(d) + 100 * max(gamma) * sd(d) * c(-1, 1), alternative = "less", j = j, alpha = alpha.up)$root, error = function(e) Inf))

    inner.ci <- cbind(inner.ci.low, inner.ci.up)
    rownames(inner.ci) <- gamma
    colnames(inner.ci) <- c("low", "up")

    outer.ci <- cbind(outer.ci.low, outer.ci.up)
    rownames(outer.ci) <- gamma
    colnames(outer.ci) <- c("low", "up")

    return(list(point.estimate = inner.ci,
                ci = outer.ci))

}

#' Compute sensitivity value
#' @param d a vector or matrix of treatment-minus-control differences (each column correponds to a hypothesis)
#' @param alpha significance level
#' @param mm test statistic, either a vector of length 3 or a matrix of three rows where each column corresponds to a U-statistic. Default is the (approximate) Wilcoxon's signed rank test.
#' @param score.method either approximate score or exact score
#' @param alternative report p-value corresponds to the maximum ("upper") or minimum ("lower") bound
#'
#' @importFrom stats qnorm
#'
#' @return sensitivity value, i.e. the kappa value such that the p-value becomes just insignificant. If \code{mm} is a matrix, then return a vector of sensitivity values corresponding to each column of \code{mm}.
#' @export
#'
#' @details The alternative direction is the the center of \code{d} is greater than 0.
#' @references Qingyuan Zhao. On sensitivity value of pair-matched observational studies. arXiv 1702.03442, \url{https://arxiv.org/abs/1702.03442}.
#'
#' @author Qingyuan Zhao
#'
#' @examples
#' d <- rnorm(100) + 1
#' gamma.star <- kappa2gamma(sen.value(d, alpha = 0.05, mm = matrix(c(2, 2, 2, 8, 5, 8), ncol = 2)))
#' gamma.star
#' sen(d, mm = c(2, 2, 2), gamma = gamma.star[1])$p.value # should equal the significance level 0.05
#'
sen.value <- function(d, alpha = 0.05, mm = c(2, 2, 2),
                      alternative = c("greater", "less", "two.sided"),
                      score.method = c("approximate", "exact")) {

    alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
    if (alternative == "less") {
        d <- -d
    } else if (alternative == "two.sided") {
        return(pmax(sen.value(d, alpha/2, mm, "greater", score.method),
                    sen.value(d, alpha/2, mm, "less", score.method)))
    }

    if (ncol(as.matrix(d)) > 1) {
        return(apply(d, 2, sen.value, alpha = alpha, mm = mm, score.method = score.method))
    }
    score.method <- match.arg(score.method, c("approximate", "exact"))
    mm <- as.matrix(mm)
    stopifnot(nrow(mm) == 3)
    r <- apply(mm, 2, function(m1) multrnks(rank(abs(d)), m1, score.method))
    T <- colSums((d > 0) * r) / colSums(r)
    c <- qnorm(1 - alpha)^2 * length(d) * colSums(r^2) / colSums(r)^2
    I <- length(d)
    kappa <- ((2 * I * T + c) - sqrt(4 * c * I * (T - T^2) + c^2)) / 2 / (I + c)
    mm.names <- paste0("(", apply(mm, 2, paste, collapse = ","), ")")
    names(kappa) <- mm.names
    return(kappa)
}

#' Power of sensitivity analysis
#'
#' @inheritParams sen.value
#' @param mu.F mean of the signed rank statistic
#' @param sigma.F standard deviation of the signed rank statistic
#' @param d empirical data used to estimate \code{mu.F} and \code{sigma.F} by jackknife
#' @param gamma target sensitivity level
#' @param alpha target significance level
#' @param I sample size
#' @param approx.method which approximation method to use?
#'
#' @return power of the sensitivity analysis, possibly a vector if \code{mm} has multiple columns.
#'
#' @details If \code{approx.method} is "fixed.alpha", then the significance level \code{alpha} is considered fixed and the corresponding quantile negligible. Otherwise we also use the \code{alpha}-quantile in the approximation formula. For more detail, see the reference.
#' @references Qingyuan Zhao. On sensitivity value of pair-matched observational studies. arXiv 1702.03442, \url{https://arxiv.org/abs/1702.03442}.
#'
#' @export
#'
#' @examples
#'
#' power.sen(d = rnorm(100) + 0.5, I = 200, gamma = 2)
#'
#' ## The following code reproduces an example of power analysis in Zhao (2017)
#' power.sen(0.76, sqrt(0.26), gamma = 2.5, I = 200)
#' power.sen(0.76, sqrt(0.26), gamma = 2.5, I = 200, approx.method = "fixed.alpha")
#'
power.sen <- function(mu.F = 1/2, sigma.F = sqrt(1/3), d = NULL, mm = c(2, 2, 2), gamma = 1, alpha = 0.05, I = 100, approx.method = c("changing.alpha", "fixed.alpha"), score.method = c("approximate", "exact")) {

    approx.method <- match.arg(approx.method, c("changing.alpha", "fixed.alpha"))
    score.method <- match.arg(score.method, c("approximate", "exact"))
    mm <- as.matrix(mm)
    stopifnot(nrow(mm) == 3)

    if (!is.null(d)) {
        mu <- function(d) {
            r <- apply(mm, 2, function(m1) multrnks(rank(abs(d)), m1, score.method))
            r <- t(t(r) / colSums(r))
            sr <- (d > 0) * r
            colSums(sr)
        }
        mu.F <- mu(d)
        mu.jackknife <- sapply(1:length(d), function(j) mu(d[-j]))
        if (ncol(mm) > 1) {
            sigma.F <- sqrt(apply(mu.jackknife, 1, var) * (length(d) - 1)^2)
        } else {
            sigma.F <- sqrt(var(mu.jackknife) * (length(d) - 1)^2)
        }
    }

    nn <- 10000
    r <- apply(mm, 2, function(m1) multrnks(1:10000, m1, score.method))
    sigma.q <- sqrt(nn * colSums(r^2) / colSums(r)^2)
    eta <- qnorm(1 - alpha)^2 * sigma.q^2 / I

    if (approx.method == "fixed.alpha") {
        mu <- mu.F - (1 / sqrt(I)) * sigma.q * qnorm(1 - alpha) * sqrt(mu.F * (1 - mu.F))
        sigma <- sigma.F / sqrt(I)
    } else {
        mu <- mu.F - ((2 * mu.F - 1) * eta + sqrt(4 * eta * mu.F * (1 - mu.F) + eta^2)) / (2 * (1 + eta))
        sigma <- (1 / sqrt(I)) * sigma.F / (1 + eta) * (1 + eta * (2 * mu.F - 1)/sqrt(4 * eta * mu.F * (1 - mu.F) + eta^2))
    }

    power <- 1 - pnorm(gamma2kappa(gamma), mu, sigma)
    if (is.null(d)) {
        return(power)
    } else {
        return(list(power = power,
                    mu.F = mu.F,
                    sigma.F = sigma.F))
    }

}

#' Bonferroni's correction with fixed \eqn{\Gamma}
#'
#' @param d a matrix of treatment-minus-control differences.
#' @param gamma sensitivity parameter (maximum odds different from a randomized experiment).
#' @inheritParams sen.value
#' @param two.sided whether a two-sided test should be used. If FALSE, test the one-sided alternative that the center of \code{d} is positive.
#'
#' @importFrom stats p.adjust
#' @return a vector of sensitivity values for each column of \code{d}
#' @export
#'
#' @details If \code{mm} is a matrix, this function computes a one-sided or two-sided p-value with each statistic (i.e. there is a p-value for every column of \code{d} and every column of $mm$), then does a Bonferroni correction over all the p-values.
#'
#' @author Qingyuan Zhao
#'
bonferroni.fg <- function(d, gamma = 1, mm = c(2, 2, 2), two.sided = TRUE) {

    I <- dim(d)[1]
    k <- dim(d)[2]
    mm <- as.matrix(mm)
    stopifnot(nrow(mm) == 3)

    p <- sapply(1:k, function(i) sen(d[, i], mm, gamma)$p.value)
    if (two.sided) {
        p <- rbind(p, sapply(1:k, function(i) sen(- d[, i], mm, gamma)$p.value))
    }
    if (k == 1) {
        p <- t(p)
    }

    if (class(p) == "matrix") {
        if (ncol(p) > 1) {
            p <- pmin(apply(p, 2, min) * nrow(p), 1)
        } else {
            p <- min(min(p) * length(p), 1)
        }
    }

    p.adjust(p, "bonferroni")

}

#' Transform sensitivity parameter in different scales
#'
#' @param kappa \eqn{\kappa = \gamma / (1 + \gamma)}
#' @export
#'
kappa2gamma <- function(kappa) {
    kappa / (1 - kappa)
}

#' @describeIn kappa2gamma Transform a sensitivity parameter from \eqn{\gamma} scale to \eqn{\kappa} scale
#'
#' @param gamma the odds of treatment of two matched units can differ at most by a factor of \code{gamma}
#' @export
#'
gamma2kappa <- function(gamma) {
    gamma / (1 + gamma)
}

#' Cross-screening
#'
#' @description Main functions that implements the cross-screening method in observational studies. \code{cross.screen} sorts the hypotheses by their sensitivity values and \code{cross.screen.fg} sorts by p-values at a fixed sensitivity \eqn{\Gamma}.
#'
#' @param d1 screen/test sample (treatment-minus-control differences), can be a matrix (rows are observations, columns are hypotheses)
#' @param d2 test/screen sample, can be a matrix
#' @param gamma sensitivity parameter (maximum odds different from a randomized experiment)
#' @param mm a vector of matrix. If matrix, adaptively choose statistic. NULL means Wilcoxon's signed rank statistic.
#' @param screen.value either "sen" (using sensitivity value) or "p" (using p-value).
#' @param alpha.screen significance level used in screening.
#' @param gamma.screen screening threshold, default is 0, meaning no screening is used.
#' @param two.sided if TRUE, automatically select the sign to test; if FALSE, test the one-sided alternative that the center of \code{d} is positive.
#' @param screen.method either keep all hypotheses significant at \code{gamma.screen} (option "threshold") or keep the least sensitive hypotheses (option "least.sensitive").
#' @param least.sensitive the number of least sensitive hypotheses to keep
#'
#' @return \code{cross.screen} returns a list
#' \describe{
#' \item{s1.kappa}{kappa values used to screen the hypotheses calculated using the first sample}
#' \item{s1.stat}{test statistics chosen using the first sample, if \code{mm} has more than 1 column}
#' \item{s1.side}{signs of alternative hypotheses chosen using the first sample}
#' \item{s1.order}{order of the hypotheses by \code{s1.kappa} if \code{s1.kappa} is above the threshold \code{gamma.screen}}
#' \item{p1}{p-values computed using the first sample at sensitivity \code{gamma}}
#' \item{s2.kappa}{kappa values used to screen the hypotheses calculated using the second sample}
#' \item{s2.stat}{test statistics chosen using the second sample, if \code{mm} has more than 1 column}
#' \item{s2.side}{signs of alternative hypotheses chosen using the second sample}
#' \item{s2.order}{order of the hypotheses by \code{s1.kappa} if \code{s1.kappa} is above the threshold \code{gamma.screen}}
#' \item{p2}{p-values computed using the second sample at sensitivity \code{gamma}}
#' \item{p}{Bonferroni adjusted p-values at sensitivity \code{gamma} computed using \code{p1} and \code{p2} (they can be directly used to control FWER)}
#' }
#' @export
#'
#' @examples
#'
#' n <- 100
#' p <- 20
#' d <- matrix(rnorm(n * p), n, p)
#' d[, 1] <- d[, 1] + 2
#' d1 <- d[1:(n/2), ]
#' d2 <- d[(n/2+1):n, ]
#' cross.screen(d1, d2,
#'              gamma = 9,
#'              gamma.screen = 1.25)$p
#'
#' ## One can run the hidden function CrossScreening:::table5(no.sims = 1)
#' ## to generate Table 5 in the paper.
#'
#' @references Qingyuan Zhao, Dylan S. Small, Paul R. Rosenbaum. Cross-screening in observational studies that test many hypotheses. arXiv preprint arXiv:1703.02078
#'
cross.screen <- function(d1,
                         d2,
                         gamma = 1,
                         mm = c(2, 2, 2),
                         screen.value = c("sen", "p"),
                         screen.method = c("threshold", "least.sensitive"),
                         alpha.screen = 0.05,
                         gamma.screen = gamma,
                         least.sensitive = 2,
                         two.sided = TRUE) {

    screen.value <- match.arg(screen.value, c("sen", "p"))
    screen.method <- match.arg(screen.method, c("threshold", "least.sensitive"))

    if (screen.value == "p") {
        cross.screen.fg(d1, d2, gamma, mm, screen.method, alpha.screen, least.sensitive, alpha.screen, gamma.screen)
    }

    k <- ncol(d1)

    ## Obtain sensitivity value for both samples and both directions of alternative
    s1.kappa.pos <- sapply(1:k, function(i) sen.value(d1[, i], alpha.screen, mm))
    s1.kappa.neg <- sapply(1:k, function(i) sen.value(- d1[, i], alpha.screen, mm))
    s2.kappa.pos <- sapply(1:k, function(i) sen.value(d2[, i], alpha.screen, mm))
    s2.kappa.neg <- sapply(1:k, function(i) sen.value(- d2[, i], alpha.screen, mm))

    if (two.sided) {
        s1.side <- 2 * (s1.kappa.pos > s1.kappa.neg) - 1 # 1 means positive alternative, -1 means negative alternative
        s1.kappa <- pmax(s1.kappa.pos, s1.kappa.neg)
        s2.side <- 2 * (s2.kappa.pos > s2.kappa.neg) - 1
        s2.kappa <- pmax(s2.kappa.pos, s2.kappa.neg)
    } else {
        s1.side <- rep(1, k)
        s1.kappa <- s1.kappa.pos
        s2.side <- rep(1, k)
        s2.kappa <- s2.kappa.pos
    }

    if (length(mm) == 3) { ## transpose everything if there is only one statistic, we want a row vector
        s1.kappa <- t(s1.kappa)
        s1.side <- t(s1.side)
        s2.kappa <- t(s2.kappa)
        s2.side <- t(s2.side)
        mm <- matrix(mm, 3, 1)
    }

    s1.stat <- apply(s1.kappa, 2, which.max)
    s1.kappa <- apply(s1.kappa, 2, max)
    s1.kappa.side <- rep(1, k)

    s2.stat <- apply(s2.kappa, 2, which.max)
    s2.kappa <- apply(s2.kappa, 2, max)
    s2.kappa.side <- rep(1, k)

    if (two.sided) {
        s1.kappa.side <- sapply(1:k, function(i) s1.side[s1.stat[i], i])
        s2.kappa.side <- sapply(1:k, function(i) s2.side[s2.stat[i], i])
    }

    ## Up till this point we have computed a sensitivity value for each hypothesis along
    ## with the corresponding statistic (if mm has more than 1 column) and side (if two.sided = TRUE)

    ## Next we order the hypotheses and remove those below the gamma.screen threshold
    if (screen.method == "threshold") {
        s1 <- order(s1.kappa, decreasing = TRUE)[1:sum(s1.kappa > gamma2kappa(gamma.screen))]
        s2 <- order(s2.kappa, decreasing = TRUE)[1:sum(s1.kappa > gamma2kappa(gamma.screen))]
    } else if (screen.method == "least.sensitive") {
        s1 <- order(s1.kappa, decreasing = TRUE)[1:least.sensitive]
        s2 <- order(s2.kappa, decreasing = TRUE)[1:least.sensitive]
    }

    ## make sure we have a hypothesis to test
    n1 <- length(s1)
    n2 <- length(s2)
    if (n1 == 0) {
        s1 <- which.max(s1.kappa)
        n1 <- 1
    }
    if (n2 == 0) {
        s2 <- which.max(s2.kappa)
        n2 <- 1
    }

    ## sensitivity analysis using the other sample
    p1 <- rep(NA, k)
    p2 <- rep(NA, k)
    if (n1 > 0) {
        p1[s1] <- sapply(s1, function(i) sen(d2[, i] * s1.kappa.side[i], mm[, s1.stat[i]], gamma = gamma)$p.value)
    }
    if (n2 > 0) {
        p2[s2] <- sapply(s2, function(i) sen(d1[, i] * s2.kappa.side[i], mm[, s2.stat[i]], gamma = gamma)$p.value)
    }

    p <- pmin(p1 * 2 * n1, p2 * 2 * n2, na.rm = TRUE)
    p <- pmin(p, 1)

    return(list(s1.kappa = s1.kappa,
                s1.stat = s1.stat,
                s1.kappa.side = s1.kappa.side,
                s1.order = s1,
                p1 = p1,
                s2.kappa = s2.kappa,
                s2.stat = s2.stat,
                s2.kappa.side = s2.kappa.side,
                s2.order = s2,
                p2 = p2,
                p = p))

}

#' @describeIn cross.screen Cross-screening with fixed \eqn{\Gamma}
#'
#' @inheritParams cross.screen
#' @return \code{cross.screen.fg} returns a list
#' \describe{
#' \item{s1.p}{p-values used to screen the hypotheses calculated using the first sample}
#' \item{s1.stat}{test statistics chosen using the first sample, if \code{mm} has more than 1 column}
#' \item{s1.side}{signs of alternative hypotheses chosen using the first sample}
#' \item{s1.order}{order of the hypotheses by \code{s1.p} if \code{s1.p} is below the threshold \code{alpha.screen}}
#' \item{p1}{p-values computed using the first sample at sensitivity \code{gamma}}
#' \item{s2.p}{p-values used to screen the hypotheses calculated using the second sample}
#' \item{s2.stat}{test statistics chosen using the second sample, if \code{mm} has more than 1 column}
#' \item{s2.side}{signs of alternative hypotheses chosen using the second sample}
#' \item{s2.order}{order of the hypotheses by \code{s2.p} if \code{s2.p} is above the threshold \code{alpha.screen}}
#' \item{p2}{p-values computed using the second sample at sensitivity \code{gamma}}
#' \item{p}{Bonferroni adjusted p-values at sensitivity \code{gamma} computed using \code{p1} and \code{p2} (they can be directly used to control FWER)}
#' }
#'
#' @export
#'
#' @examples
#'
#' ## The following code generates Table 1 in the paper.
#'
#' require(CrossScreening)
#' data(nhanes.fish)
#' data(nhanes.fish.match)
#'
#' data <- nhanes.fish
#' match <- nhanes.fish.match
#'
#' outcomes <- grep("^o\\.", names(data))
#' log2diff <- function(y1, y2) {
#'     if (min(c(y1, y2)) == 0) {
#'         y1 <- y1 + 1
#'         y2 <- y2 + 1
#'     }
#'     log2(y1) - log2(y2)
#' }
#' d <- sapply(outcomes, function(j) log2diff(data[match$treated, j], data[match$control, j]))
#' set.seed(11)
#' split <- sample(1:nrow(d), nrow(d) / 2, replace = FALSE)
#' d1 <- d[split, ]
#' d2 <- d[-split, ]
#'
#' mm <- matrix(c(2, 2, 2, 8, 5, 8), ncol = 2)
#' data.frame(outcome = names(data)[outcomes],
#'            p.value =
#'                cross.screen(d1, d2,
#'                             gamma = 9,
#'                             screen.value = "p",
#'                             screen.method = "least.sensitive",
#'                             mm = mm)$p)
#'
#'
#' @author Qingyuan Zhao
#'
cross.screen.fg <- function(d1, d2,
                            gamma = 1,
                            mm = c(2, 2, 2),
                            screen.method = c("threshold", "least.sensitive"),
                            alpha.screen = 0.05,
                            gamma.screen = gamma,
                            least.sensitive = 2,
                            two.sided = TRUE) {

    k <- ncol(d1)
    if (!is.null(mm)) {
        mm <- as.matrix(mm)
    }
    screen.method <- match.arg(screen.method, c("threshold", "least.sensitive"))

    ## Obtain sensitivity value for both samples and both directions of alternative
    p1.screen.pos <- sapply(1:k, function(i) sen(d1[, i], mm, gamma = gamma.screen)$p.value)
    p1.screen.neg <- sapply(1:k, function(i) sen(- d1[, i], mm, gamma = gamma.screen)$p.value)
    p2.screen.pos <- sapply(1:k, function(i) sen(d2[, i], mm, gamma = gamma.screen)$p.value)
    p2.screen.neg <- sapply(1:k, function(i) sen(- d2[, i], mm, gamma = gamma.screen)$p.value)

    if (two.sided) {
        p1.side <- 2 * (p1.screen.pos < p1.screen.neg) - 1 # 1 means positive alternative, -1 means negative alternative
        p1.screen <- pmin(p1.screen.pos, p1.screen.neg)
        p2.side <- 2 * (p2.screen.pos < p2.screen.neg) - 1
        p2.screen <- pmin(p2.screen.pos, p2.screen.neg)
    } else {
        p1.side <- rep(1, k)
        p1.screen <- p1.screen.pos
        p2.side <- rep(1, k)
        p2.screen <- p2.screen.pos
    }

    if (is.null(mm) || (ncol(mm) == 1)) {
        p1.screen <- t(p1.screen)
        p1.side <- t(p1.side)
        p2.screen <- t(p2.screen)
        p2.side <- t(p2.side)
    }

    p1.screen.statistic <- apply(p1.screen, 2, which.min)
    p1.screen <- apply(p1.screen, 2, min)
    p1.screen.side <- rep(1, k)

    p2.screen.statistic <- apply(p2.screen, 2, which.min)
    p2.screen <- apply(p2.screen, 2, min)
    p2.screen.side <- rep(1, k)

    if (two.sided) {
        p1.screen.side <- sapply(1:k, function(i) p1.side[p1.screen.statistic[i], i])
        p2.screen.side <- sapply(1:k, function(i) p2.side[p2.screen.statistic[i], i])
    }

    if (screen.method == "threshold") {
        s1 <- order(p1.screen)[1:sum(p1.screen <= alpha.screen)]
        s2 <- order(p2.screen)[1:sum(p2.screen <= alpha.screen)]
    } else if (screen.method == "least.sensitive") {
        s1 <- order(p1.screen)[1:least.sensitive]
        s2 <- order(p2.screen)[1:least.sensitive]
    }

    ## make sure we have a hypothesis to test
    n1 <- length(s1)
    n2 <- length(s2)
    if (n1 == 0) {
        s1 <- which.min(p1.screen)
        n1 <- 1
    }
    if (n2 == 0) {
        s2 <- which.min(p2.screen)
        n2 <- 1
    }

    ## p reported
    p1 <- rep(NA, k)
    p2 <- rep(NA, k)
    p1[s1] <- sapply(s1, function(i) sen(d2[, i] * p1.screen.side[i],
                                         mm[, p1.screen.statistic[i]],
                                         gamma)$p.value)
    p2[s2] <- sapply(s2, function(i) sen(d1[, i] * p2.screen.side[i],
                                         mm[, p1.screen.statistic[i]],
                                         gamma)$p.value)
    p <- pmin(p1 * 2 * n1,
              p2 * 2 * n2, 1)

    return(list(s1.p = p1.screen,
                s1.stat = p1.screen.statistic,
                s1.side = p1.screen.side,
                s1.order = s1,
                p1 = p1,
                s2.p = p2.screen,
                s2.stat = p2.screen.statistic,
                s2.side = p2.screen.side,
                s2.order = s2,
                p2 = p2,
                p = p))

}

#' Fallback procedure for multiple testing
#'
#' @param p a vector of p-values
#' @param alpha significance level
#' @param spread the way to spread \code{alpha}, either a vector of the same length as \code{p} or a single number to indicate equal spread in the first \code{spread} hypotheses.
#'
#' @return the rejected hypotheses (TRUE means reject, FALSE means accept)
#' @export
#'
#' @author Qingyuan Zhao
#' @references Brian L. Wiens. A fixed sequence Bonferroni procedure for testing multiple endpoints. Pharmaceutical Statistics, 2(3), 211---215, 2003.
#'
fallback.test <- function(p, alpha = 0.05, spread = 1) {

    if (length(p) == 0) {
        return(NULL)
    }

    if (sum(!is.na(p)) == 0) {
        return(rep(NA, length(p)))
    }

    if (length(spread) == 1) {
        spread <- min(spread, sum(!is.na(p)))
        alpha.seq <- c(rep(alpha / spread, spread), rep(0, length(p) - spread))
    } else {
        alpha.seq <- spread
    }

    alpha.current <- 0
    reject <- rep(0, length(p))
    for (i in 1:length(p)) {
        reject[i] <- as.numeric(p[i] <= alpha.current + alpha.seq[i])
        if (!is.na(reject[i]) & (reject[i] == 1)) {
            alpha.current <- alpha.current + alpha.seq[i]
        } else {
            alpha.current <- 0
        }
    }

    return(reject == 1)
}

#' Recycling procedure for multiple testing
#'
#' @inheritParams fallback.test
#'
#' @return rejected hypotheses
#'
#' @details WARNING: only supports recycle the first two tests.
#' @export
#'
#' @author Qingyuan Zhao
#'
recycle.test <- function(p, alpha = 0.05) {

    ## handle missing values
    if (length(p) == 0) {
        return(NULL)
    }

    if (sum(!is.na(p)) == 0) {
        return(rep(NA, length(p)))
    }

    if (sum(!is.na(p)) == 1){
        return(as.numeric(p <= alpha))
    }

    reject <- rep(0, length(p))

    alpha.current <- alpha / 2
    if (p[1] <= alpha.current) {
        reject[1] <- 1
        alpha.current <- alpha
    } else {
        alpha.current <- alpha / 2
    }

    if (p[2] <= alpha.current) {
        reject[2] <- 1
        if (reject[1] == 0) { ## recycle
            if (p[1] <= alpha) {
                reject[1] <- 1
                alpha.current <- alpha
            } else {
                alpha.current <- 0
            }
        } else {
            alpha.current <- alpha.current
        }
    } else {
        alpha.current <- 0
    }

    if (length(p) >= 3) {
        reject[3:length(p)] <- fallback.test(p[3:length(p)], alpha.current, 1)
    }

    return(reject == 1)

}
