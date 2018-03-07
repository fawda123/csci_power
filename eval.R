library(tidyverse)
library(readxl)
library(nlme)

fls <- list.files('R', full.names = T)
sapply(fls, source)

# raw <- read_excel("ignore/csci_raw.xlsx")
# save(raw, file = 'data/raw.RData', compress = 'xz')

data(raw)

# format raw data for power analysis, CSCI
topow <- raw %>% 
  select(MasterID, FieldReplicate, SampleYear, CSCI) %>% 
  filter(SampleYear >= 2008) %>% 
  filter(!is.na(CSCI)) %>% 
  mutate(SampleYear = as.character(SampleYear))

# tmp <- topow %>% 
#   group_by(MasterID) %>% 
#   mutate(n = n())
# 
# ggplot(tmp, aes(x = n)) + 
#   geom_bar()

## 
# input for function

# indicator name
indicator <- 'CSCI'

# mean of indicator
ind.mean <- topow$CSCI %>% 
  mean(na.rm = T)

# power to detect trend of x percent
trend <- 2 

# nsites - number of sites visited each year
# a matrix where number of rows is number of rows is number of indicators
# number of columns is number of years
# each cell is number of sites to visit in a year
# this is 30 sites visited every year
# skipped years should be zero
nyr <- 34 # up to 2050 (from 2016)
nsites1 <- matrix(rep(60, nyr), nrow = 1, ncol = nyr)
nsites2 <- matrix(rep(c(50, 0), times = 10), nrow = 1, ncol = 20)
nsites3 <- matrix(rep(c(30, 0), times = nyr), nrow = 1, ncol = nyr)
nsites4 <- matrix(rep(c(1, 1), times = nyr), nrow = 1, ncol = nyr)
nsites <- rbind(nsites1, nsites2, nsites3, nsites4)

# number of repeats per year (one)
nrepeats <- 1

##
# variance components of sites

mod <- lme(CSCI ~ SampleYear, random = ~ 1 | MasterID, data = topow)
tot.var <- var(topow$CSCI)

year.var <- mod$sigma^2
index.var <- var(resid(mod))
site.var <- tot.var - year.var - index.var

power.fcn(indicator = 'CSCI', 
          ind.mean = ind.mean, 
          trend = trend,
          nsites = nsites2, 
          nrepeats = 1, 
          site.var = site.var, 
          year.var = year.var, 
          siteyear.var = 0, 
          index.var = index.var, 
          site.rho = 1, 
          year.rho = 0, 
          alfa = 0.05, 
          plot.ind = T)
