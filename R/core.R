#' Sensivity analysis with signed score test
#'
#' @param y a vector of treatment-minus-control differences
#' @param mm a vector (m, munder, mover) that indicates the U-statistic. NULL means Wilcoxon's signed rank test.
#' @param gamma sensitivity parameter >= 1.
#' @param tail report p-value corresponds to the maximum ("upper") or minimum ("lower") bound
#'
#' @importFrom sensitivitymw multrnks
#' @importFrom stats pnorm
#'
#' @return A list
#' \describe{
#' \item{pval}{p-value}
#' \item{pval2}{two sided p-value}
#' \item{T}{test statistic}
#' \item{E}{Mean of the test statistic under sensivity gamma}
#' \item{V}{Variance of the test statistic under sensitivity gamma}
#' \item{devc}{Effect size of T compared to E and V}
#' }
#'
#' @export
#' @author Paul Rosenbaum, Qingyuan Zhao
#'
sen <- function(y, mm = NULL, gamma = 1, tail = c("upper", "lower")) {

    tail <- match.arg(tail, c("upper", "lower"))
    if (tail == "lower") {
        y <- - y
    }

    ## compute score based on the rank of abs(y)
    if (is.null(mm)) {
        qs <- rank(abs(y))
    } else {
        mm <- as.matrix(mm)
        stopifnot(nrow(mm) == 3)
        if (ncol(mm) == 1) {
            qs <- multrnks(rank(abs(y)), m1 = mm[2], m2 = mm[3], m = mm[1])
        } else {
            sen.result <- apply(mm, 2, function(mmm) sen(y, mmm, gamma))
            return(data.frame(do.call("rbind", lapply(sen.result, unlist))))
        }
    }
    qs <- qs * (abs(y) > 0)
    sg <- as.numeric(y > 0) # sign
    Ts <- sum(qs * sg) # test statistic
    kappa <- gamma / (1 + gamma)
    ETGamma1 <- sum(qs)/2 # expectation of the statistic when Gamma = 1
    VT <- sum(qs*qs)*kappa*(1-kappa) # variance of the bounding statistic
    ET <- sum(qs) * kappa
    dev <- (Ts - ET) / sqrt(VT)
    pval <- 1 - pnorm(dev)
    kappahat <- Ts / sum(qs)
    gammahat <- max(kappahat / (1 - kappahat), (1 - kappahat) / kappahat)

    pval2 <- min(2 * pval,1)
    devc <- max(0, dev)

    return(list(pval=pval,
                pval2=pval2,
                gammahat=gammahat,
                kappahat=kappahat,
                T=Ts,
                E=ET,
                V=VT,
                devc=devc))
}


#' Compute sensitivity value
#' @param d a vector of treatment-minus-control differences
#' @param alpha significance level
#' @param mm test statistic, either a vector of length 3 or a matrix of three rows where each column corresponds to a U-statistic. Default is the (approximate) Wilcoxon's signed rank test.
#'
#' @importFrom sensitivitymw multrnks
#' @importFrom stats qnorm
#'
#' @return sensitivity value, i.e. the kappa value such that the p-value becomes just insignificant. If \code{mm} is a matrix, then return a vector of sensitivity values corresponding to each column of \code{mm}.
#' @export
#'
#' @details The alternative direction is the the center of \code{d} is greater than 0.
#'
#' @references Qingyuan Zhao. On sensitivity value of pair-matched observational studies. arXiv 1702.03442, \url{https://arxiv.org/abs/1702.03442}.
#'
#' @author Qingyuan Zhao
#'
get.kappa <- function(d, alpha = 0.05, mm = c(2, 2, 2)) {
    mm <- as.matrix(mm)
    stopifnot(nrow(mm) == 3)
    r <- apply(mm, 2, function(m1) multrnks(rank(abs(d)), m1[2], m1[3], m1[1]))
    T <- colSums((d > 0) * r) / colSums(r)
    c <- qnorm(1 - alpha)^2 * length(d) * colSums(r^2) / colSums(r)^2
    I <- length(d)
    kappa <- ((2 * I * T + c) - sqrt(4 * c * I * (T - T^2) + c^2)) / 2 / (I + c)
    return(kappa)
}

#' Bonferroni's correction with fixed \eqn{\Gamma}
#'
#' @param d a matrix of treatment-minus-control differences.
#' @param gamma sensitivity parameter (maximum odds different from a randomized experiment).
#' @inheritParams get.kappa
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

    p <- sapply(1:k, function(i) sen(d[, i], mm, gamma)$pval)
    if (two.sided) {
        p <- rbind(p, sapply(1:k, function(i) sen(- d[, i], mm, gamma)$pval))
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

#' Transform a sensitivity parameter from kappa scale to gamma scale
#'
#' @keywords internal
#'
kappa.to.gamma <- function(kappa) {
    kappa / (1 - kappa)
}

#' Transform a sensitivity parameter from gamma scale to kappa scale
#'
#' @keywords internal
#'
gamma.to.kappa <- function(gamma) {
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
#' @param alpha.screen significance level used in screening.
#' @param gamma.screen screening threshold, default is 0, meaning no screening is used.
#' @param two.sided if TRUE, automatically select the sign to test; if FALSE, test the one-sided alternative that the center of \code{d} is positive.
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
                         alpha.screen = 0.05,
                         gamma.screen = 0,
                         two.sided = TRUE) {

    k <- ncol(d1)

    ## Obtain sensitivity value for both samples and both directions of alternative
    s1.kappa.pos <- sapply(1:k, function(i) get.kappa(d1[, i], alpha.screen, mm))
    s1.kappa.neg <- sapply(1:k, function(i) get.kappa(- d1[, i], alpha.screen, mm))
    s2.kappa.pos <- sapply(1:k, function(i) get.kappa(d2[, i], alpha.screen, mm))
    s2.kappa.neg <- sapply(1:k, function(i) get.kappa(- d2[, i], alpha.screen, mm))

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
    s1 <- order(s1.kappa, decreasing = TRUE)[1:sum(s1.kappa > gamma.to.kappa(gamma.screen))]
    s2 <- order(s2.kappa, decreasing = TRUE)[1:sum(s1.kappa > gamma.to.kappa(gamma.screen))]

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
        p1[s1] <- sapply(s1, function(i) sen(d2[, i] * s1.kappa.side[i], mm[, s1.stat[i]], gamma = gamma)$pval)
    }
    if (n2 > 0) {
        p2[s2] <- sapply(s2, function(i) sen(d1[, i] * s2.kappa.side[i], mm[, s2.stat[i]], gamma = gamma)$pval)
    }

    p <- pmin(p1 * 2 * n1,
              p2 * 2 * n2, 1)

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
#' @param screen.method either keep all hypotheses significant at \code{gamma.screen} (option "threshold") or keep the least sensitive hypotheses (option "least_sensitive").
#' @param alpha.least.sensitive the number of least sensitive hypotheses to keep
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
#'                cross.screen.fg(d1, d2,
#'                                gamma = 9,
#'                                screen.method = "least_sensitive",
#'                                mm = mm)$p)
#'
#'
#' @author Qingyuan Zhao
#'
cross.screen.fg <- function(d1, d2,
                            gamma = 1,
                            screen.method = c("threshold", "least_sensitive"),
                            alpha.screen = 0.05,
                            alpha.least.sensitive = 2,
                            gamma.screen = gamma,
                            mm = c(2, 2, 2),
                            two.sided = TRUE) {

    k <- ncol(d1)
    if (!is.null(mm)) {
        mm <- as.matrix(mm)
    }
    screen.method <- match.arg(screen.method, c("threshold", "least_sensitive"))

    ## Obtain sensitivity value for both samples and both directions of alternative
    p1.screen.pos <- sapply(1:k, function(i) sen(d1[, i], mm, gamma = gamma.screen)$pval)
    p1.screen.neg <- sapply(1:k, function(i) sen(- d1[, i], mm, gamma = gamma.screen)$pval)
    p2.screen.pos <- sapply(1:k, function(i) sen(d2[, i], mm, gamma = gamma.screen)$pval)
    p2.screen.neg <- sapply(1:k, function(i) sen(- d2[, i], mm, gamma = gamma.screen)$pval)

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
    } else if (screen.method == "least_sensitive") {
        s1 <- order(p1.screen)[1:alpha.least.sensitive]
        s2 <- order(p2.screen)[1:alpha.least.sensitive]
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
                                         gamma)$pval)
    p2[s2] <- sapply(s2, function(i) sen(d1[, i] * p2.screen.side[i],
                                         mm[, p1.screen.statistic[i]],
                                         gamma)$pval)
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
#' @return the rejected hypotheses
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

    return(reject)
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

    reject

}
