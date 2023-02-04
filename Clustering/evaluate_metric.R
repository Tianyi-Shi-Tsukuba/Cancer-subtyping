library(survival)

survfile <- paste("GLIO_Survival.txt", sep="")
surv <- read.table(survfile, header = T)
labels=clust$clustering
survresult<-survdiff(Surv(Survival, Death)~labels, data=surv)
options(digtis = 5)
p.val <- 1 - pchisq(survresult$chisq, length(survresult$n) - 1)
print(p.val)
print(-log10(p.val))
