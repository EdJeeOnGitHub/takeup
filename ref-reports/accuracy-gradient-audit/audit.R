#!/usr/bin/env Rscript
# Descriptive reporting audit, not a structural refit or a causal distance effect.
out <- "ref-reports/accuracy-gradient-audit"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source("R/structural/asymmetric-observability-data.R")
e <- new.env()
load("build/structural-workspace/main-core-input.RData", e)
s <- e$stan_data
# This checks respondent-level counts against the saved full structural input.
checked <- main_core_prepare_peer_response_data(s)
a <- s$analysis_data
cl <- unique(a[, c("cluster.id", "cluster.dist.to.pot", "dispersed_community")])
stopifnot(!anyDuplicated(cl$cluster.id))
p <- readRDS("data/clean-data/clean-endline-know-table-data-long.rds")
p <- p[as.character(p$know.table.type) == "table.A" &
         as.character(p$sms.treatment) == "sms.control" &
         p$cluster.id %in% cl$cluster.id, ]
p$recognized_flag <- as.character(p$recognized) == "yes"
nr <- tapply(p$recognized_flag, p$KEY.individ, sum)
p <- p[p$KEY.individ %in% names(nr)[nr > 0], ]
stopifnot(nrow(p) == s$num_beliefs_obs * s$know_table_A_sample_size)
p$distance <- cl$cluster.dist.to.pot[match(p$cluster.id, cl$cluster.id)] / 1000
p$dispersed <- cl$dispersed_community[match(p$cluster.id, cl$cluster.id)]
p$arm <- factor(as.character(p$assigned.treatment),
                levels = c("control", "calendar", "ink", "bracelet"))
p$truth <- ifelse(p$actual.other.dewormed.any.1, "participant", "nonparticipant")
p$definite <- as.character(p$dewormed) %in% c("yes", "no")
p$correct <- (as.character(p$dewormed) == "yes") == p$actual.other.dewormed.any.1
slopes <- predictions <- tests <- counts <- list()
for (sample in c("all_communities", "exclude_dispersed")) {
  q <- if (sample == "all_communities") p else p[!p$dispersed, ]
  counts[[length(counts)+1L]] <- data.frame(sample,
    communities = length(unique(q$cluster.id)), respondents = length(unique(q$KEY.individ)),
    peer_rows = nrow(q), linked_rows = sum(!is.na(q$truth)),
    recognized_linked = sum(q$recognized_flag & !is.na(q$truth)),
    definite_linked = sum(q$recognized_flag & !is.na(q$truth) & q$definite))
  q <- q[q$recognized_flag & !is.na(q$truth), ]
  for (stage in c("definite", "accuracy")) {
    d <- if (stage == "accuracy") q[q$definite, ] else q
    d$y <- if (stage == "accuracy") as.integer(d$correct) else as.integer(d$definite)
    for (truth in c("nonparticipant", "participant")) {
      z <- d[d$truth == truth, ]
      fit <- glm(y ~ 0 + arm + arm:distance, family = binomial(), data = z)
      stopifnot(fit$converged, all(is.finite(coef(fit))))
      V <- sandwich::vcovCL(fit, cluster = z$cluster.id, type = "HC1")
      b <- coef(fit)
      ix <- grep(":distance$", names(b))
      G <- length(unique(z$cluster.id))
      crit <- qt(.975, G-1)
      se <- sqrt(diag(V))[ix]
      slopes[[length(slopes)+1L]] <- data.frame(sample, stage, truth,
        arm = levels(z$arm), n = as.integer(table(z$arm)), clusters = G,
        log_odds_per_km = b[ix], se, lower = b[ix]-crit*se, upper = b[ix]+crit*se)
      for (restriction in c("all_slopes_zero", "equal_slopes_across_arms")) {
        L <- diag(length(b))[ix, , drop=FALSE]
        if (restriction == "equal_slopes_across_arms") L <- L[-1,,drop=FALSE] - L[rep(1,3),,drop=FALSE]
        v <- L %*% b
        W <- as.numeric(t(v) %*% solve(L %*% V %*% t(L), v))
        df <- nrow(L)
        tests[[length(tests)+1L]] <- data.frame(sample, stage, truth, restriction,
          F = W/df, df1 = df, df2 = G-1, p = pf(W/df,df,G-1,lower.tail=FALSE))
      }
      for (arm in levels(z$arm)) {
        nd <- data.frame(arm=factor(rep(arm,2),levels=levels(z$arm)),distance=c(.5,2.5))
        X <- model.matrix(delete.response(terms(fit)), nd)
        pr <- plogis(drop(X %*% b))
        grad <- pr*(1-pr)*X
        gd <- grad[2,]-grad[1,]
        ds <- sqrt(drop(t(gd)%*%V%*%gd))
        predictions[[length(predictions)+1L]] <- data.frame(sample, stage, truth, arm,
          p_500m=pr[1],p_2500m=pr[2],difference=pr[2]-pr[1],
          lower=pr[2]-pr[1]-crit*ds,upper=pr[2]-pr[1]+crit*ds)
      }
    }
  }
}
write.csv(do.call(rbind,counts),file.path(out,"sample-counts.csv"),row.names=FALSE)
write.csv(do.call(rbind,slopes),file.path(out,"slopes.csv"),row.names=FALSE)
write.csv(do.call(rbind,predictions),file.path(out,"predictions.csv"),row.names=FALSE)
write.csv(do.call(rbind,tests),file.path(out,"joint-tests.csv"),row.names=FALSE)
plot_data <- subset(do.call(rbind,predictions),
  stage=="accuracy" & sample=="all_communities")
pdf(file.path(out,"accuracy-gradients.pdf"),width=9,height=4.5)
par(mfrow=c(1,2),mar=c(5,6,3,1))
for (truth in c("nonparticipant","participant")) {
  z <- plot_data[plot_data$truth==truth,]
  plot(100*z$difference,4:1,xlim=c(-35,65),ylim=c(.5,4.5),yaxt="n",
    xlab="Change from 0.5 to 2.5 km (percentage points)",ylab="",
    main=paste("True",truth),pch=19)
  axis(2,at=4:1,labels=z$arm,las=1)
  segments(100*z$lower,4:1,100*z$upper,4:1)
  abline(v=0,lty=2,col="gray50")
}
dev.off()
writeLines(trimws(capture.output(sessionInfo()), which="right"), file.path(out,"session-info.txt"))
print(do.call(rbind,counts))
print(subset(do.call(rbind,tests),stage=="accuracy"),row.names=FALSE)
print(subset(do.call(rbind,predictions),stage=="accuracy" & sample=="all_communities"),row.names=FALSE)
