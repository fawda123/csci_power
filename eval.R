library(tidyverse)
library(readxl)
library(nlme)
library(doParallel)
library(foreach)

# raw <- read_excel("ignore/csci_raw.xlsx")
# save(raw, file = 'data/raw.RData', compress = 'xz')

data(raw)

# format raw data for power analysis, CSCI
topow <- raw %>% 
  select(MasterID, FieldReplicate, SampleYear, CSCI) %>% 
  filter(SampleYear >= 2008 & SampleYear <= 2015) %>% 
  filter(!is.na(CSCI)) %>% 
  mutate(yrno = SampleYear - min(SampleYear))

# tmp <- topow %>%
#   group_by(MasterID) %>%
#   mutate(n = n())
# 
# ggplot(tmp, aes(x = n)) +
#   geom_bar()

## 
# input for function

# input mixed model
mod <- lme(CSCI ~ yrno, random = ~ 1 | MasterID, data = topow)

scns <- crossing(
  nsite = c(15, 30, 60),
  dec = seq(0.01, 0.05, length = 5),
  L = 10:35,
  s.freq = 1:2
  )

library(doParallel)
ncores <- detectCores() - 1  

# setup parallel backend
cl<-makeCluster(ncores)
registerDoParallel(cl)
strt<-Sys.time()

# process all stations
res <- foreach(i = 1:nrow(scns), .packages = 'nlme') %dopar% {
  
  sink('log.txt')
  cat(i, 'of', nrow(scns), '\n')
  print(Sys.time()-strt)
  sink()
  
  source("R/funcs.R")
  
  dec <- scns[i, ][['dec']]
  nsite <- scns[i,][['nsite']]
  L <- scns[i, ][['L']]
  s.freq <- scns[i, ][['s.freq']]
  
  power.fun(mod, dec = dec, s.freq = s.freq, nsite = nsite, L = L, boot.num = 2000)

}
save(res, file = 'data/res.RData', compress = 'xz')  
pows <- lapply(res, function(x) x$p.boot) %>% unlist

toplo <- scns %>% 
  mutate(
    pow = 100 * pows,
    dec = factor(dec, levels = seq(0.01, 0.05, length = 5), labels = paste(seq(1, 5, length = 5), '%')), 
    L = 2015 + L, 
    nsite = paste(nsite, 'sites'), 
    s.freq = factor(s.freq, levels = c(1, 2), labels = c('annual', 'biennial'))
    )

p <- ggplot(toplo, aes(x = L, y = pow, group = dec, colour = factor(dec))) +
  # geom_line() +
  geom_point(colour = NA, guide = F) +
  geom_smooth(se = F) +
  facet_grid(s.freq ~ nsite) +
  theme_bw(base_family = 'serif') +
  scale_y_continuous('Power') +
  scale_colour_manual('Trend', values = RColorBrewer::brewer.pal(9, 'Reds')[4:9]) +
  theme(
    strip.background = element_blank(), 
    legend.position = 'top', 
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

png('power.png', height = 6, width = 8, units = 'in', res = 400)
p
dev.off()