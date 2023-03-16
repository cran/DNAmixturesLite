## Analysing trace MC18, using USCaucasian allele-frequencies.
data(MC18, USCaucasian, SGMplusDyes)
## The prosecution hypothesis could be that K1, K2, and K3 have
## contributed to the trace. Setting detection threshold to 50rfu.
## For layout as an EPG we follow the layout of SGMplus.
mixHp <- DNAmixture(list(MC18), k = 3, K = c("K1", "K2", "K3"), C = list(50),
                  database = USCaucasian, dyes = list(SGMplusDyes))
plot(mixHp, epg = TRUE)
## The yellow is a bit difficult to see; we can
## change the colors representing the dyes
plot(mixHp, epg = TRUE, dyecol = list(c("blue", "green", "black")))

## Fit by maximum likelihood using p as the initial value for optimisation
p <- mixpar(rho = list(30), eta = list(34), xi = list(0.08),
            phi = list(c(K1 = 0.71, K3 = 0.1, K2 = 0.19)))
mlHp <- mixML(mixHp, pars = p)
mlHp$mle

## Find the estimated covariance matrix of the MLE
V.Hp <- varEst(mixHp, mlHp$mle, npars = list(rho=1,eta=1,xi=1,phi=1))
V.Hp$cov ## using (rho, eta)
## The summary is a table containing the MLE and their standard errors
summary(V.Hp)

## Assess the fit of the distribution of peak heights
qqpeak(mixHp, pars = mlHp$mle, dist = "conditional")

## Assess the prediction of allele presence
pr <- prequential.score(mixHp, pars = mlHp$mle)
plot(pr)

## Compare observed peak heights to peak heights simulated under the model
sims <- rPeakHeight(mixHp, mlHp$mle, nsim = 100, dist = "conditional")
oldpar <- par(mfrow = c(2,5), mar = c(2,2,2,0))
boxplot(mixHp, sims)
par(oldpar)

\donttest{
data(MC18, USCaucasian, SGMplusDyes)

## A defence hypothesis could be that one unknown person has
## contributed along with K1 and K3.
mixHd <- DNAmixture(list(MC18), k = 3, K = c("K1", "K3"), C = list(50),
                  database = USCaucasian, dyes = list(SGMplusDyes))

## Fit by ML
mlHd <- mixML(mixHd, pars = mixpar(rho = list(30), eta = list(34), xi = list(0.08),
                       phi = list(c(K1 = 0.71, K3 = 0.1, U1 = 0.19))))
mlHd$mle

## Find the estimated covariance matrix of the MLE
V <- varEst(mixHd, mlHd$mle, npars = list(rho=1,eta=1,xi=1,phi=1))
summary(V)

## Assess the peak height distribution
qq <- qqpeak(mixHd, pars = mlHd$mle, dist = "conditional")

## Assess the prediction of allele presence
pr <- prequential.score(mixHd, pars = mlHd$mle)
plot(pr, col = pr$trace)

## Compare observed peak heights to peak heights simulated under the model
sims <- rPeakHeight(mixHd, mlHd$mle, nsim = 100, dist = "conditional")
par(mfrow = c(2,5), mar = c(2,2,2,0))
boxplot(mixHd, sims)


## The log10 of likelihood-ratio can be found as
(mlHp$lik - mlHd$lik)/log(10)


## We can find the most probable genotypes for U1 under the hypothesis
## Hd : K1, K2, and U1.

## Include the peak heights
setPeakInfo(mixHd, mlHd$mle)
## Jointly best for all unknown contributors
mp <- map.genotypes(mixHd, pmin = 0.01, type = "seen")
## See the genotypes rather than allelecounts
print(summary(mp), marker = "D2S1338")

## Include only information on presence and absence of alleles
setPeakInfo(mixHd, mlHd$mle, presence.only = TRUE)
## Jointly best for all unknown contributors
mp2 <- map.genotypes(mixHd, pmin = 0.01, type = "seen")
## See the genotypes rather than allelecounts
print(summary(mp2), marker = "D2S1338")
}
