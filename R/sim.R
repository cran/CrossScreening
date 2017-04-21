#' Generate Table 5 of the paper
#'
#' @param no.sims number of simulations
#'
#' @import plyr
#' @import parallel
#' @import tables
#'
#' @details Be careful, this starts a parallel program that may use all your CPU!
#' @import stats
#'
#' @keywords internal
#'
table5 <- function(no.sims = 1) {

    warning("Be careful, if no.sims > 1, this starts a parallel program that may use all your CPU!")

    ## simulation settings
    no.outcomes.vec <- c(100)
    no.subjects.vec <- c(100, 250, 500)
    alpha.screen.vec <- c(0.05)
    error.vec <- c("normal")
    signal.vec <- c(1, 2, 3)
    zeta.vec <- c(0.8)
    effect.vec <- c(0.5)
    mm.list <- list(c(2, 2, 2),
                    c(8, 5, 8),
                    matrix(c(8, 5, 8, 8, 6, 7, 8, 7, 8), ncol = 3))
    gamma.vec <- c(2)
    settings <- expand.grid(no.outcomes = no.outcomes.vec,
                            no.subjects = no.subjects.vec,
                            alpha.screen = alpha.screen.vec,
                            error = error.vec,
                            signal = signal.vec,
                            zeta = zeta.vec,
                            effect = effect.vec,
                            gamma = gamma.vec,
                            gamma.screen = gamma.vec,
                            mm = c(1, 2, 3),
                            two.sided = c(TRUE))

    no.outcomes.vec <- c(500)
    no.subjects.vec <- c(500)
    alpha.screen.vec <- c(0.05)
    error.vec <- c("normal")
    signal.vec <- c(1, 2, 3)
    zeta.vec <- c(0.8)
    effect.vec <- c(0.5)
    mm.list <- list(c(2, 2, 2),
                    c(8, 5, 8),
                    matrix(c(8, 5, 8, 8, 6, 7, 8, 7, 8), ncol = 3))
    gamma.vec <- c(2)
    settings <- rbind(settings,
                      expand.grid(no.outcomes = no.outcomes.vec,
                                  no.subjects = no.subjects.vec,
                                  alpha.screen = alpha.screen.vec,
                                  error = error.vec,
                                  signal = signal.vec,
                                  zeta = zeta.vec,
                                  effect = effect.vec,
                                  gamma = gamma.vec,
                                  gamma.screen = gamma.vec,
                                  mm = c(1, 2, 3),
                                  two.sided = c(TRUE)))



    one.simulation <- function(setting) {

        d <- switch(as.character(setting$error),
                    normal = rnorm(setting$no.outcomes*setting$no.subjects),
                    t4 = rt(setting$no.outcomes*setting$no.subjects, 4),
                    t3 = rt(setting$no.outcomes*setting$no.subjects, 3))
        d <- matrix(d, ncol = setting$no.outcomes)
        ## effect <- switch(as.character(setting$error),
        ##                  normal = 1,
        ##                  t4 = sqrt(2),
        ##                  t3 = sqrt(3))
        ## effect <- effect * setting$effect.size
        if (setting$signal == 1) {
            effect <- switch(as.character(setting$error),
                             normal = c(0.5),
                             t4 = c(0.7))
        } else if (setting$signal == 2) {
            effect <- switch(as.character(setting$error),
                             normal = c(0.5, 0.5),
                             t4 = c(0.7, 0.7))
        } else if (setting$signal == 3) {
            effect <- switch(as.character(setting$error),
                             normal = c(0.6, 0.4),
                             t4 = c(0.8, 0.6))
        }
        d[, 1:length(effect)] <- t(t(d[, 1:length(effect)]) + effect)

        mm.this <- mm.list[[setting$mm]]

        ## Bonferroni
        p.bonf <- bonferroni.fg(d, setting$gamma, mm.this, setting$two.sided)
        reject.bonf <- which(p.bonf < 0.05)

        ## Cross screening
        s.half <- setting$no.subjects / 2
        p.cross <- cross.screen(d[1:s.half, ],
                               d[(s.half + 1):setting$no.subjects, ],
                               gamma = setting$gamma,
                               alpha.screen = setting$alpha.screen,
                               gamma.screen = 0,
                               mm = mm.this,
                               two.sided = setting$two.sided)
        ## Order the hypotheses being tested
        o1 <- p.cross$s1.order
        o2 <- p.cross$s2.order
        ## Cross screening: Bonferroni
        oo1 <- which(p.cross$s1.kappa > gamma2kappa(setting$gamma.screen))
        oo2 <- which(p.cross$s2.kappa > gamma2kappa(setting$gamma.screen))
        n1 <- length(oo1)
        n2 <- length(oo2)
        if (n1 > 0) {
            reject1 <- oo1[which(p.cross$p1[oo1] * n1 < 0.05 / 2)]
        } else {
            reject1 <- integer(0)
        }
        if (n2 > 0) {
            reject2 <- oo2[which(p.cross$p2[oo2] * n2 < 0.05 / 2)]
        } else {
            reject2 <- integer(0)
        }
        reject.cross.bonf <- union(reject1, reject2)
        ## Cross screening: fixed sequence
        reject1 <- o1[which(fallback.test(p.cross$p1[o1], alpha = 0.05 / 2) == 1)]
        reject2 <- o2[which(fallback.test(p.cross$p2[o2], alpha = 0.05 / 2) == 1)]
        reject.cross.fixed <- union(reject1, reject2)
        ## Cross screening: fallback (alpha/2, alpha/2)
        reject1 <- o1[which(fallback.test(p.cross$p1[o1], alpha = 0.05 / 2, spread = 2) == 1)]
        reject2 <- o2[which(fallback.test(p.cross$p2[o2], alpha = 0.05 / 2, spread = 2) == 1)]
        reject.cross.fallback <- union(reject1, reject2)
        ## Cross screening: recycle
        reject1 <- o1[which(recycle.test(p.cross$p1[o1], alpha = 0.05/2) == 1)]
        reject2 <- o2[which(recycle.test(p.cross$p2[o2], alpha = 0.05/2) == 1)]
        reject.cross.recycle <- union(reject1, reject2)

        ## Single splitting
        s <- round(setting$no.subjects * (1 - setting$zeta))
        p.single <- cross.screen(d[1:s, ],
                                d[(s + 1):setting$no.subjects, ],
                                gamma = setting$gamma,
                                alpha.screen = setting$alpha.screen,
                                gamma.screen = setting$gamma.screen,
                                mm = mm.this,
                                two.sided = setting$two.sided)
        ## Single splitting: fixed sequence
        o1 <- p.single$s1.order
        reject.single.fixed <- o1[which(fallback.test(p.single$p1[o1], alpha = 0.05) == 1)]

        reject <- rbind(c(1, 2) %in% reject.bonf,
                        c(1, 2) %in% reject.cross.bonf,
                        c(1, 2) %in% reject.cross.fixed,
                        c(1, 2) %in% reject.cross.fallback,
                        c(1, 2) %in% reject.cross.recycle,
                        c(1, 2) %in% reject.single.fixed)
        colnames(reject) <- c("H1", "H2")

        output <- cbind(setting, reject, row.names = NULL)
        output$method <- c("Bonferroni", "Cross Bonf", "Cross fixed", "Cross fallback", "Cross recycle", "Single fixed")

        output

    }

    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores)

    clusterExport(cl, ls(pos = globalenv()))
    output <- parLapply(cl,
                        1:no.sims,
                        function(sim) { do.call("rbind", lapply(1:nrow(settings), function(i) {one.simulation(settings[i, ])})) })

    stopCluster(cl)

    output <- do.call(rbind, output)

    no.outcomes <- no.subjects <- alpha.screen <- error <- signal <- zeta <- effect <- gamma <- mm <- two.sided <- method <- H1 <- H2 <- NULL # This is just to avoid "Undefined global functions or variables" NOTE from "R CMD check"
    output.methods <- setdiff(unique(output$method), "Cross Bonf")
    output.summary <- ddply(subset(output, method %in% output.methods),
                            .(no.outcomes, no.subjects, alpha.screen, error, signal, zeta, effect, gamma, mm, two.sided, method),
                            summarise, H1.power = mean(H1), H2.power = mean(H2), H12.power = mean(H1 & H2))
    output.summary$no.outcomes <- factor(output.summary$no.outcomes)
    output.summary$no.subjects <- factor(output.summary$no.subjects)
    output.summary$signal <- factor(output.summary$signal)
    output.summary$statistic <- factor(output.summary$mm)
    output.summary$method <- factor(output.summary$method)
    output.summary$gamma <- factor(output.summary$gamma)
    output.summary$two.sided <- factor(output.summary$two.sided)
    output.summary$H1.power <- round(output.summary$H1.power * 100, 1)
    output.summary$H2.power <- round(output.summary$H2.power * 100, 1)
    output.summary$H12.power <- round(output.summary$H12.power * 100, 1)
    levels(output.summary$signal) <- c("(0.5,0)", "(0.5,0.5)", "(0.6,0.4)")
    levels(output.summary$signal) <- c("(0.7,0)", "(0.7,0.7)", "(0.8,0.6)")
    levels(output.summary$statistic) <- c("Wilcoxon", "(8,5,8)", "Adapative")
    output.summary$H2.power[output.summary$signal == "(0.5,0)"] <- NA
    output.summary$H12.power[output.summary$signal == "(0.5,0)"] <- NA

    df.100 <- subset(output.summary, no.outcomes == 100)
    df.100 <- droplevels(df.100)
    df.500 <- subset(output.summary, no.outcomes == 500)
    df.500 <- droplevels(df.500)

    table.100 <- tabular(no.outcomes * no.subjects * signal * statistic ~ method * (H1.power * Format(digits = 2) + H2.power * Format(digits = 2) + H12.power * Format(digits = 2)) * Heading() * identity, data = subset(df.100))
    print(table.100)

    table.500 <- tabular(no.outcomes * no.subjects * signal * statistic ~ method * (H1.power * Format(digits = 3) + H2.power * Format(digits = 2) + H12.power * Format(digits = 2)) * Heading() * identity, data = subset(df.500))
    print(table.500)

}
