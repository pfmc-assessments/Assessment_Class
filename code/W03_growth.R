## ---- include=FALSE--------------------------------------------------------------------------------------------------------
library(ggplot2)
library(magrittr)
theme_set(theme_classic(base_size = 16))




## ---- eval=FALSE-----------------------------------------------------------------------------------------------------------
## library(nwfscSurvey)
## rex.bio <- PullBio.fn(Name = "rex sole",
##                       SurveyName = "NWFSC.Combo",
##                       SaveFile = FALSE)




## ----growth_scatter, fig.show='hide'---------------------------------------------------------------------------------------
ggplot(rex.bio) +
  geom_point(aes(x = Age, y = Length_cm, col = Sex), 
             alpha = 0.25) 




## --------------------------------------------------------------------------------------------------------------------------
unsexed.index <- which(rex.bio$Sex=='U' & !is.na(rex.bio$Age))
set.seed(20398)
males <- sample(unsexed.index, 
                size = floor(length(unsexed.index)/2), 
                replace = FALSE)
rex.bio2 <- rex.bio
rex.bio2$Sex[males] <- 'M'
rex.bio2$Sex[setdiff(unsexed.index, males)] <- 'F'


## --------------------------------------------------------------------------------------------------------------------------
rex.bio2$Sex <- factor(rex.bio2$Sex)

vbgf.nls <- nls(Length_cm ~ linf[Sex]*
                  (1-exp(-k[Sex]*(Age-t0[Sex]))), 
                data = rex.bio2, 
                start = list(linf = rep(40,2), 
                             k = rep(0.2,2), 
                             t0 = rep(1,2)))


## --------------------------------------------------------------------------------------------------------------------------
summary(vbgf.nls)


## --------------------------------------------------------------------------------------------------------------------------
summary(vbgf.nls, cor=TRUE)$correlation %>%
  round(2)


## --------------------------------------------------------------------------------------------------------------------------
vbgf.nls2 <- nls(Length_cm ~ la1[Sex] + (la2[Sex] - la1[Sex]) * 
                   (1-exp(-k[Sex]*(Age-1))) / 
                   (1-exp(-k[Sex]*19)), 
                 data = rex.bio2, 
                 start = list(la1 = rep(12,2), 
                              la2 = rep(35,2), 
                              k = rep(0.2,2)))


## --------------------------------------------------------------------------------------------------------------------------
summary(vbgf.nls2)


## --------------------------------------------------------------------------------------------------------------------------
summary(vbgf.nls2, cor=TRUE)$correlation %>%
  round(2)


## ---- echo=FALSE-----------------------------------------------------------------------------------------------------------
plot(rex.bio2$Age[!is.na(rex.bio2$Age)], resid(vbgf.nls2),
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
dat.m <- dplyr::filter(rex.bio2, Sex == 'M', !is.na(Age))
dat.f <- dplyr::filter(rex.bio2, Sex == 'F', !is.na(Age))

pars.init <- log(c(la1.m = 13, la1.f = 13, 
                   la2.m = 35, la2.f = 38, 
                   k.m = .15, k.f = .17, cv = .1))

vbgf.optim <- optim(pars.init, vbgf.loglik, 
                    dat.m = dat.m, dat.f = dat.f, 
                    a1 = 1, a2 = 20)



## ---- results='hold'-------------------------------------------------------------------------------------------------------
exp(vbgf.optim$par)[1:4]
exp(vbgf.optim$par)[5:7]

