# mod.in from lme
# L is number of years
# dec is proportion decline from beginnign to end of L
# boot.num is number of random mods to eval
# s.freq is sampling frequency, 1 for every year, 2 for biannual, etc.
# alpha is alpha val for two-tailed test
# nsite is how many places you look

# assumes the likelihood of observing a trend is the ability to show a trend across all sampled sites
power.fun<-function(mod.in,L=20,dec=0.5,boot.num=1000,s.freq=1,alpha=0.05, nsite = 30){
  
  # get starting value, intercept from model (make sure model starts at year zero)
  No<-as.numeric(fixef(mod.in))[1]
  
  # indicator name
  indic <- formula(mod.in$terms)
  indic <- as.character(indic[2])
  
  ##
  # variance components
  # total variation is sum of year, site, and index variation
  # year variation is sigma from model
  # index variation is residual variance from model
  # site variation is remainder
  
  # total variation
  tot.var <- var(mod.in$data[[indic]])
  
  # annual variation
  year.var <- mod$sigma^2
  index.var <- var(resid(mod))
  site.var <- tot.var - year.var - index.var
  
  #a is decline per year based on total in dec from 1:L
  a<--1*dec*No/L
  
  # get simulations for one to many sites
  temp.data<-TS.sim.rawreg(L, No, a, index.var, site.var, boot.num * nsite)
  temp.data <- matrix(temp.data, nrow = L * nsite)
  
  #gets model coefficients relating random model values to time
  #if significant, indicates sig decrease over period of obs
  reg.test<-function(dat, yrs, subs, sitein) {
    tomod <- data.frame(dat, yrs, sitein)[subs, ]
    tomod <- data.frame(tomod, subs)
    # mod <- lme(dat~ yrs, random = ~ 1 | sitein, data = tomod)
    # sums <- summary(mod)$tTable
    # sums[2,5]
    mod <- lm(dat~ yrs + sitein, data = tomod)
    summary(mod)$coefficients[2, 4]
  }

  #power calcs
  subs<-seq(1,L * nsite, by=s.freq)
  yrs <- rep(seq(1, L), times = nsite)
  sites <- factor(rep(1:nsite, each = L))
  p.out<-apply(temp.data,2,reg.test, yrs = yrs, subs = subs, sitein = sites)
  
  #number of models with sig change divided total num of models
  #i.e., number of models in which a change was detected - power
  p.reg.reject<-length(which(p.out<alpha))/boot.num
  
  #output
  out.ls<-list(p.reg.reject,No,L,a,unique(yrs),nsite)
  names(out.ls)<-c('p.boot','start.val','obs.period','chg.per.yr',
                   'samp.yrs','nsite')
  out.ls
  
}

##
#dependent functions
#this gets predicted values based on uncertainty over period of obs.
TS.sim.rawreg<-function(L, No, a, index.var, site.var, boot.num)  {
  yr<-1:L
  N<-No + yr*a - a
  rand.var<-rnorm(L*boot.num,0,index.var) + rnorm(L*boot.num,0,site.var)
  rand.var<-matrix(rand.var,nrow=L,ncol=boot.num)
  apply(rand.var,2,function(x) x+N)
}
