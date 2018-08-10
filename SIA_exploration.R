## explore whisker sia data ##
rm(list = ls())
library(tidyverse)
library(lubridate)
library(ggpubr)
options(dplyr.width = Inf) #enables head() to display all coloums

# load data
# adults analysed were regrowths
# pup whiskers were plucked at the end of fieldwork
d0 <- read_csv("./data/lnfs_all_whisker_bulk_sia_correctedsequence.csv")
d1 <- d0 %>% 
  filter(adult == TRUE) %>% 
  select(seal, num.from.base, c13, n15, sample.len, root, root.length, original.length, plucked, start.date, end.date, glycine) %>% 
  mutate(start.date = as.Date(strptime(start.date, format = '%d/%m/%y')), end.date = as.Date(strptime(end.date, format = '%d/%m/%y'))) %>% 
  arrange(seal, num.from.base)

# quick look at isotope values ----
d1 %>% 
  ggplot(aes(y = n15, x = seal, group = seal, colour = glycine)) + 
  geom_boxplot()

d1 %>% 
  ggplot(aes(y = c13, x = seal, group = seal, colour = glycine)) + 
  geom_boxplot()

d1 %>% 
  group_by(glycine) %>% 
  summarise(maxN = max(n15, na.rm = T)) # max d15N is 16 for KI shelf area (agrees with Lowther et al. 2013)

d1 %>% 
  ggplot(aes(x = c13, y = n15, colour = glycine)) + 
  geom_point() + 
  facet_wrap(~seal)


# find cutoff n15 value for glycine seals ----
g <- d1 %>% filter(glycine == TRUE)
g %>% ggplot(aes(y = n15, x = num.from.base)) + geom_point(size = 0.5) + geom_line() + facet_wrap(~seal) + 
  geom_hline(aes(yintercept = 17.6))

# label glycine segments
d1 <- d1 %>% mutate(has.glycine = ifelse(n15 > 17.6, TRUE, FALSE))
d1$has.glycine[is.na(d1$has.glycine)] <- FALSE

# estimate dates of whisker segments ----
# for adult regrowths
interp.regrowth.dates <- function(df){
  # cols needed: 
  # start.date = date whisker was first cut, end.date = date regrowth was collected, 
  # original length = length of regrowth, plucked = if whisker was plucked with root or cut, 
  # root.length = length of root (under the skin), len.from.base = length from base of regrowth in mm
  # has.glycine = logical if segment represents glycine values
  x <- df
  x <- x %>% 
    arrange(num.from.base) %>% 
    mutate(num.from.tip = rev(num.from.base))
  
  # non-glycine whiskers ----
  if(x$glycine[1] == FALSE){
    # assign grouping for the regrowth section 
    # get row number corresponding to start.date by accounting for root length from the tip. 
    l <- x %>% 
      arrange(num.from.tip) %>% 
      mutate(len.from.tip = cumsum(sample.len)) %>% .$len.from.tip
    rn <- which.min(abs(x$root.length[1] - l)) + 1
    
    x1 <- x %>% 
      arrange(num.from.tip) %>% 
      mutate(rg = 0)
    x1$rg[rn:nrow(x1)] <- seq(1,nrow(x1)-rn+1,1)
    
    # calculate growth.rate and resolution ---
    if(x$plucked[1] == TRUE){
      # plucked non-glycine whiskers
      x <- x1 %>% mutate(res = round(difftime(end.date,start.date, units = 'days')/(original.length - root.length), 2),
                         growth.rate = (original.length - root.length)/as.integer(difftime(end.date, start.date, units = 'days')))
    } else {
      # cut non-glycine whiskers
      x <- x1 %>% mutate(res = round(difftime(end.date,start.date, units = 'days')/original.length, 2),
                         growth.rate = original.length/as.integer(difftime(end.date, start.date, units = 'days')))
    }
  } else { 
    # plucked glycine whiskers ----
    # filter out only regrowth portion based on glycine values. 
    x <- x %>% 
      arrange(num.from.tip) %>% 
      mutate(rg = cumsum(has.glycine))
    
    x1 <- x %>% filter(rg > 0)
    
    # calculate growth.rate and temporal resolution
    x1 <- x1 %>% mutate(res = round(difftime(end.date, start.date, units = 'days')/sum(sample.len), 2), 
                        growth.rate = sum(sample.len)/as.integer(difftime(end.date, start.date, units = 'days')))
    res <- x1$res[1]
    growth.rate <- x1$growth.rate[1]
    x <- x %>% mutate(res = res, growth.rate = growth.rate)
  }
  
  # estimate dates 
  # starting from end.date (base) ---
  # x <- x %>%
  #   arrange(num.from.base) %>% 
  #   # calculate cumulative length from base
  #   mutate(len.from.base = cumsum(sample.len)) %>% 
  #   # estimate dates started from the end.date (newest info)
  #   arrange(len.from.base) %>% 
  #   mutate(date.est = end.date - (lag(len.from.base) * res))
  # x$date.est[1] <- x$end.date[1]
  
  # starting from start.date ---
  x1 <- x %>% 
    filter(rg > 0) %>% 
    arrange(num.from.tip) %>% 
    mutate(len.from.start = cumsum(sample.len)) %>%
    # estimate dates started from the end.date (newest info)
    mutate(date.est = start.date + (lag(len.from.start) * res))
  
  # then estimates dates for segments before start.date
  x2 <- full_join(x, x1) %>% filter(is.na(date.est))
  if(nrow(x2) > 1) {
    x2 <- x2 %>% 
    arrange(num.from.base) %>% 
    mutate(len.from.start = cumsum(sample.len)) %>%
    # estimate dates before the start.date
    mutate(date.est = start.date - (lag(len.from.start) * res)) %>% 
    mutate(len.from.start = len.from.start  * -1 + first(len.from.start)) %>% 
    # filter(!is.na(date.est)) %>% 
    arrange(num.from.tip)
    x2$len.from.start[is.na(x2$date.est)] <- x2$sample.len[is.na(x2$date.est)]
    x2$date.est[is.na(x2$date.est)] <- x2$start.date[1]
  x <- bind_rows(x2, x1) %>% filter(!is.na(date.est))
  } else {
    x1$date.est[1] <- x1$start.date[1]
    x <- x1
  }
  x %>% arrange(num.from.base)
}

d2 <- d1 %>% 
  group_by(seal) %>% 
  nest() %>% 
  mutate(data = map(data, interp.regrowth.dates)) %>% 
  unnest()

# look at growth rates 
d2 %>% 
  group_by(seal) %>% 
  summarise(growth.rate = first(growth.rate), glycine = first(glycine)) # 450b dates are likely the correct ones

# remove seal 450, keep seal , rename 450b to 450
d2 <- d2 %>% filter(!seal == '450') 
d2$seal[d2$seal == '450b'] <- '450'


# Plots ----
# create sequence from tip (oldest info) 
d2 <- d2 %>% 
  arrange(seal, num.from.base) %>% 
  group_by(seal) %>% 
  mutate(num.from.tip = rev(num.from.base))

# modify estimated dates for plotting purposes
d2 <- d2 %>% 
  mutate(date = as.Date(format(date.est, format = '2017-%m-%d'))) %>% 
  mutate(date = ifelse(date > ymd(20171101), 
                       format(date, format = '2016-%m-%d'), 
                       format(date, format = '2017-%m-%d'))) %>% 
  mutate(date = as.Date(date)) %>% 
  mutate(start.date = as.Date(format(start.date, format = '2017-%m-%d')), end.date = as.Date(format(end.date, format = '2017-%m-%d')))
  

# plot c13 and n15 time series
pdf(file = './plots/SI_plots.pdf', width = 13, height = 11)
pal <- wesanderson::wes_palette('BottleRocket2', 3, type = 'continuous' )
d2 %>% 
  ggplot(aes(x = c13, y = n15, colour = num.from.tip)) +
  geom_point() +
  facet_wrap(glycine~seal) + 
  # scale_color_distiller(type = 'div')
  scale_colour_gradient2(low = pal[1], mid = pal[2], high = pal[3], midpoint = 13) +
  theme_bw()

d2 %>% 
  ggplot(aes(x = date, y = c13)) + 
  geom_point(size = 1, aes(colour = has.glycine)) + 
  geom_line(size = 0.5) +
  geom_vline(aes(xintercept = start.date)) +
  geom_vline(aes(xintercept = end.date)) +
  labs(x = 'date') +
  scale_x_date(date_breaks = '1 month', date_labels = '%b') +
  theme_bw() + 
  facet_wrap(~seal)

d2 %>% 
  filter(has.glycine == FALSE) %>% 
  ggplot(aes(x = date, y = n15)) + 
  geom_point(size = 1) + 
  geom_line(size = 0.5) +
  geom_vline(aes(xintercept = start.date)) +
  geom_vline(aes(xintercept = end.date)) +
  labs(x = 'date') +
  scale_x_date(date_breaks = '1 month', date_labels = '%b') +
  theme_bw() + 
  facet_wrap(~seal)

d2 %>% 
  filter(glycine == TRUE) %>% 
  ggplot(aes(x = date, y = n15)) + 
  geom_point(size = 1, aes(colour = has.glycine)) + 
  geom_line(size = 0.5) +
  geom_vline(aes(xintercept = start.date)) +
  geom_vline(aes(xintercept = end.date)) +
    geom_hline(aes(yintercept = 17.6), linetype = 'dashed') +
  labs(x = 'date') +
  scale_x_date(date_breaks = '1 month', date_labels = '%b') +
  theme_bw() + 
  facet_wrap(~seal) + 
  scale_colour_manual(values = wesanderson::wes_palette('Rushmore1')[c(3,5)])
dev.off()

# outliers -----
# detect outliers
d1 %>% ggplot(aes(x = seal, y = c13)) +
  geom_boxplot() 
d1 %>% ggplot(aes(x = seal, y = n15)) +
  geom_boxplot() # 326 suspiciously high n15 value. 
d1 %>% ggplot(aes(x = c13, y = n15)) + 
  geom_point() + 
  facet_wrap(~seal)

# remove outliers
x <- d1$c13[d1$seal == '103']
boxplot.stats(x)$out
outlier_values <- boxplot.stats(x)$out  # outlier values.
boxplot(x, boxwex=0.1)
mtext(paste("Outliers: ", paste(outlier_values, collapse=", ")), cex=0.6)
d1$c13[d1$seal == '103'][which(d1$c13[d1$seal == '103'] %in% outlier_values )] <- NA

# investigate delay in glycine absorption in whisker -----------------------------------------------------------
g <- d1 %>% filter(glycine == TRUE) %>% mutate(seal = paste0(seal, 'g'))
g1 <- d1 %>% filter(glycine == TRUE) %>% mutate(glycine = FALSE)
g2 <- bind_rows(g, g1)

# label glycine segments
g2 <- g2 %>% mutate(has.glycine = ifelse(n15 > 17.6, TRUE, FALSE))
g2$has.glycine[is.na(g2$has.glycine)] <- FALSE

# estimate dates of whisker segments ----
g3 <- g2 %>% 
  group_by(seal) %>% 
  nest() %>% 
  mutate(data = map(data, interp.regrowth.dates)) %>% 
  unnest()

# look at growth rates 
g3 %>% 
  group_by(seal, glycine) %>% 
  summarise(growth.rate = first(growth.rate))

# create sequence from tip (oldest info) 
g3 <- g3 %>% 
  arrange(seal, num.from.tip)

# modify estimated dates for plotting purposes
g3 <- g3 %>% 
  mutate(date = as.Date(format(date.est, format = '2017-%m-%d'))) %>% 
  mutate(date = ifelse(date > ymd(20171101), 
                       format(date, format = '2016-%m-%d'), 
                       format(date, format = '2017-%m-%d'))) %>% 
  mutate(date = as.Date(date))


# plot c13 and n15 time series
g3 %>% 
  ggplot(aes(x = date, y = c13)) + 
  geom_point(size = 1, aes(colour = has.glycine)) + 
  geom_line(size = 0.5) + 
  geom_vline(aes(xintercept = start.date)) +
  geom_vline(aes(xintercept = end.date)) +
  labs(x = 'date') +
  scale_x_date(date_breaks = '1 month', date_labels = '%b') +
  theme_bw() + 
  facet_wrap(~seal) 

d2 %>% 
  filter(has.glycine == FALSE) %>% 
  ggplot(aes(x = date, y = n15)) + 
  geom_point(size = 1) + 
  geom_line(size = 0.5) +
  labs(x = 'date') +
  scale_x_date(date_breaks = '1 month', date_labels = '%b') +
  theme_bw() + 
  facet_wrap(~seal)

