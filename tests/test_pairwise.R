source("R/pairwise_tests.R")
set.seed(741)
z <- data.frame(group = factor(rep(LETTERS[1:4], each = 35)),
  time = rexp(140, rep(c(.2, .3, .7, 1), each = 35)) + .01,
  event = rbinom(140, 1, .8), age = rnorm(140))
pw <- pairwise_logrank(z$time, z$event, z$group)
stopifnot(nrow(pw) == 6, all(pw$p_holm >= pw$p),
  isTRUE(all.equal(pw$p_holm, p.adjust(pw$p, "holm"))))
zz <- z[z$group %in% c("A", "B"), ]; zz$group <- droplevels(zz$group)
direct <- survival::survdiff(survival::Surv(time, event) ~ group, data = zz)
stopifnot(abs(pw$statistic[1] - direct$chisq) < 1e-10)
x <- round(z$time, 1)
dn <- pairwise_numeric(x, z$group)
reference <- rstatix::dunn_test(data.frame(x, group = z$group), x ~ group)
stopifnot(all(abs(dn$p - reference$p) < 1e-10), all(abs(dn$p_holm - reference$p.adj) < 1e-10))
fit <- survival::coxph(survival::Surv(time, event) ~ group + age, data = z, x = TRUE, model = TRUE)
cp <- pairwise_cox(fit, z, "group")
for (i in seq_len(nrow(cp))) {
  zz <- z; zz$group <- relevel(zz$group, ref = cp$group2[i])
  ff <- survival::coxph(survival::Surv(time, event) ~ group + age, data = zz)
  ss <- summary(ff); term <- paste0("group", cp$group1[i])
  stopifnot(abs(cp$HR[i] - exp(coef(ff)[term])) < 1e-8,
    abs(cp$p[i] - ss$coefficients[term, "Pr(>|z|)"]) < 1e-8)
}
zero <- pairwise_logrank(z$time, rep(0, nrow(z)), z$group)
stopifnot(nrow(zero) == 6, all(is.na(zero$p)), all(zero$status == "No events"))
stopifnot(nrow(pairwise_numeric(1:5, rep("one", 5))) == 0)
cat("Pairwise tests: log-rank, Dunn ties, Holm families, Cox contrasts, edge cases PASSED\n")
