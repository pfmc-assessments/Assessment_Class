## ---- include=FALSE--------------------------------------------------------------------------------------------------------
library(ggplot2)
library(magrittr)
theme_set(theme_classic(base_size = 16))




## ---- eval=FALSE-----------------------------------------------------------------------------------------------------------
## library(nwfscSurvey)
## canary.bio <- PullBio.fn(Name = "canary rockfish",
##                          SurveyName = "NWFSC.Combo",
##                          SaveFile = FALSE)




## ----growth_scatter, fig.show='hide'---------------------------------------------------------------------------------------
ggplot(canary.bio) +
  geom_point(aes(x = Age, y = Length_cm, col = Sex), 
             alpha = 0.25) 




## --------------------------------------------------------------------------------------------------------------------------
unsexed.index <- which(canary.bio$Sex=='U' & 
                         !is.na(canary.bio$Age))
set.seed(20398)
males <- sample(unsexed.index, 
                size = floor(length(unsexed.index)/2), 
                replace = FALSE)
canary.bio2 <- canary.bio
canary.bio2$Sex[males] <- 'M'
canary.bio2$Sex[setdiff(unsexed.index, males)] <- 'F'


## --------------------------------------------------------------------------------------------------------------------------
canary.bio2$Sex <- factor(canary.bio2$Sex)

vbgf.nls <- nls(Length_cm ~ linf[Sex]*
                  (1-exp(-k[Sex]*(Age-t0[Sex]))), 
                data = canary.bio2, 
                start = list(linf = rep(55,2), 
                             k = rep(0.2,2), 
                             t0 = rep(0,2)))


## --------------------------------------------------------------------------------------------------------------------------
summary(vbgf.nls)


## --------------------------------------------------------------------------------------------------------------------------
summary(vbgf.nls, cor=TRUE)$correlation %>%
  round(2)


## --------------------------------------------------------------------------------------------------------------------------
vbgf.nls2 <- nls(Length_cm ~ la1[Sex] + (la2[Sex] - la1[Sex]) * 
                   (1-exp(-k[Sex]*(Age-1))) / 
                   (1-exp(-k[Sex]*24)), 
                 data = canary.bio2, 
                 start = list(la1 = rep(12,2), 
                              la2 = rep(50,2), 
                              k = rep(0.2,2)))


## --------------------------------------------------------------------------------------------------------------------------
summary(vbgf.nls2)


## --------------------------------------------------------------------------------------------------------------------------
summary(vbgf.nls2, cor=TRUE)$correlation %>%
  round(2)


## ---- echo=FALSE-----------------------------------------------------------------------------------------------------------
plot(canary.bio2$Age[!is.na(canary.bio2$Age)], resid(vbgf.nls2),
     xlab = 'Age')
abline(h=0)


## --------------------------------------------------------------------------------------------------------------------------
vbgf.loglik <- function(log.pars, dat.m, dat.f, a1, a2) {
  pars <- exp(log.pars)
  
  l.pred.m <- pars['la1.m'] + (pars['la2.m'] - pars['la1.m']) * 
    (1-exp(-pars['k.m']*(dat.m$Age - a1))) / 
    (1-exp(-pars['k.m']*(a2-a1)))
  
  l.pred.f <- pars['la1.f'] + (pars['la2.f'] - pars['la1.f']) * 
    (1-exp(-pars['k.f']*(dat.f$Age - a1))) / 
    (1-exp(-pars['k.f']*(a2-a1)))
  
  nll <- -dlnorm(x = c(dat.m$Length_cm, dat.f$Length_cm), 
                 meanlog = log(c(l.pred.m, l.pred.f)) -
                   pars['cv']^2/2, 
                 sdlog = pars['cv'], log = TRUE) %>%
    sum()
  return(nll)
}



## --------------------------------------------------------------------------------------------------------------------------
dat.m <- dplyr::filter(canary.bio2, Sex == 'M', !is.na(Age))
dat.f <- dplyr::filter(canary.bio2, Sex == 'F', !is.na(Age))

pars.init <- log(c(la1.m = 12, la1.f = 12, 
                   la2.m = 50, la2.f = 55, 
                   k.m = .17, k.f = .14, cv = .1))

vbgf.optim <- optim(pars.init, vbgf.loglik, 
                    dat.m = dat.m, dat.f = dat.f, 
                    a1 = 1, a2 = 25)



## ---- results='hold'-------------------------------------------------------------------------------------------------------
exp(vbgf.optim$par)[1:4]
exp(vbgf.optim$par)[5:7]


## ---- echo=FALSE-----------------------------------------------------------------------------------------------------------
nls.f <- coef(vbgf.nls2)[c(1,3,5)]
nls.m <- coef(vbgf.nls2)[c(2,4,6)]
optim.f <- exp(vbgf.optim$par)[c(2,4,6)]
optim.m <- exp(vbgf.optim$par)[c(1,3,5)]

ggplot(canary.bio2) +
  geom_point(aes(x = Age, y = Length_cm, col = Sex), alpha = 0.25) +
  geom_line(aes(x = Age, y = nls.f[1] + (nls.f[2] - nls.f[1]) *
                  (1-exp(-nls.f[3]*(Age-1))) / 
                  (1-exp(-nls.f[3]*24)), linetype='1'), col = 'blue') +
  geom_line(aes(x = Age, y = nls.m[1] + (nls.m[2] - nls.m[1]) *
                  (1-exp(-nls.m[3]*(Age-1))) / 
                  (1-exp(-nls.m[3]*24)), linetype='1'), col = 'red') +
  geom_line(aes(x = Age, y = optim.f[1] + (optim.f[2] - optim.f[1]) *
                  (1-exp(-optim.f[3]*(Age-1))) / 
                  (1-exp(-optim.f[3]*24)), linetype='2'), col = 'blue') +
  geom_line(aes(x = Age, y = optim.m[1] + (optim.m[2] - optim.m[1]) *
                  (1-exp(-optim.m[3]*(Age-1))) / 
                  (1-exp(-optim.m[3]*24)), linetype='2'), col = 'red') +
  scale_color_manual(values = c('F' = 'blue', 'M' = 'red')) +
  # scale_linetype_identity(name = "Model fit",
  #                         breaks = c("solid", "dotted"),
  #                         labels = c("normal", "lognormal"),
  #                         guide = "legend") +
  NULL

