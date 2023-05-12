## explore whisker sia data ##
rm(list = ls())
library(tidyverse)
library(lubridate)
library(ggpubr)
library(viridisLite)
library(mapdata)
library(viridis)
options(dplyr.width = Inf) #enables head() to display all coloums

## Load data
# adults analysed were regrowths
# pup whiskers were plucked at the end of fieldwork
d0 <- read_csv("./data/lnfs_all_whisker_bulk_sia_correctedsequence.csv")
d1 <- d0 %>% 
  filter(adult == TRUE) %>% 
  select(seal, num.from.base, c13, n15, sample.len, root, root.length, original.length, plucked, start.date, end.date, glycine) %>% 
  mutate(start.date = as.Date(strptime(start.date, format = '%d/%m/%y')), end.date = as.Date(strptime(end.date, format = '%d/%m/%y'))) %>% 
  arrange(seal, num.from.base)

# remove seal 450, keep seal , rename 450b to 450
d1 <- d1 %>% filter(!seal == '450') 
d1$seal[d1$seal == '450b'] <- '450'

# quick look at isotope values ----
d1 %>% filter(n15 <= 17.6) %>% 
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


# plot carbon/nitrogen in sequence 
# pup data 
p <- d0 %>% filter(adult == FALSE) %>% 
  select(seal, num.from.base, c13, n15, sample.len, root, root.length, original.length, plucked, start.date, end.date, glycine, growth.rate, res) %>% 
  mutate(start.date = as.Date(strptime(start.date, format = '%d/%m/%y')), end.date = as.Date(strptime(end.date, format = '%d/%m/%y'))) %>% 
  arrange(seal, num.from.base) %>% 
  filter(c13 > -18) %>% 
  mutate(seal = ifelse(seal == '104', '104_77', '103_71')) %>% 
  mutate(has.glycine = FALSE)

# join mum and pup
mp <- bind_rows(d1, p)
mp %>%  ggplot(aes(num.from.base, c13)) + 
  geom_line() + 
  geom_point(size = 1, aes(colour = has.glycine)) + 
  facet_wrap(~seal) + 
  scale_colour_brewer(type = 'qual', palette = 'Set1')

mp %>% 
  filter(has.glycine == FALSE) %>% 
  ggplot(aes(num.from.base, n15)) + 
  geom_line() + 
  geom_point(size = 1, aes(colour = has.glycine)) + 
  facet_wrap(~seal) + 
  scale_colour_brewer(type = 'qual', palette = 'Set1')


# Glycine seals only ----------------------------------------------------------
g1 <- g %>% mutate(has.glycine = ifelse(n15 > 17.6, TRUE, FALSE))

## is the whisker length pre-deployment equal to the root length? ##
tmp <- g1 %>% group_by(seal) %>% 
  filter(has.glycine == TRUE) %>% 
  mutate(gly.start = max(num.from.base)) %>% 
  select(seal, gly.start)
tmp <- tmp[!duplicated(tmp),]

g1 <- left_join(g1, tmp)

a <- g1 %>% filter(num.from.base > gly.start) %>% 
  group_by(seal) %>% 
  summarise(predeploy.len = sum(sample.len), 
            root.length = first(root.length), 
            diff = predeploy.len - root.length, 
            ratio = predeploy.len/root.length)

a %>% filter(ratio <= 1) %>% 
  summarise(mean.ratio = mean(ratio)) # pre-deployment root length average 0.8 of current root length. 

## estimate growth rate from glycine seals 
a <- g1 %>% group_by(seal) %>% 
  group_by(seal) %>% 
  mutate(days.since.injection = as.numeric(difftime(end.date, start.date, units = 'days')), 
         len.from.base = cumsum(sample.len)) %>% 
  filter(has.glycine == T) %>% 
  mutate(max.num.from.base = max(num.from.base)) %>% 
  filter(num.from.base == max.num.from.base) %>% 
  mutate(growth.rate.yearly = (365/days.since.injection) * len.from.base,
         growth.rate.daily = growth.rate.yearly/365,
         time.res = 1/growth.rate.daily) %>% 
  select(seal, growth.rate.daily, time.res, max.num.from.base)

g2 <- left_join(g1, a) %>% filter(num.from.base <= max.num.from.base)
g2 <- g2 %>% group_by(seal) %>% 
  arrange(seal, desc(num.from.base)) %>% 
  mutate(len.from.start = cumsum(sample.len), date.est = start.date + time.res * lag(len.from.start))

g2 %>% ggplot(aes(date.est, c13)) + geom_point() + geom_line() + facet_wrap(~seal)


# estimate dates of whisker segments ----------------------------------------
# for adult regrowths
interp.regrowth.dates <- function(df, ratio = 1, method = "auto"){
  # cols needed: 
  # start.date = date whisker was first cut, 
  # end.date = date regrowth was collected, 
  # original length = length of regrowth, 
  # plucked = if whisker was plucked with root or cut, 
  # root.length = length of root (under the skin), 
  # len.from.base = length from base of regrowth in mm
  # has.glycine = logical if segment represents glycine values
  
  # "auto" method: if there's glycine. use glycine to estimate growth.rate. 
  # If there's no glycine use manual regrowth rate accounting for plucked/cut whiskers
  x <- df
  x <- x %>% 
    arrange(num.from.base) %>% 
    mutate(num.from.tip = rev(num.from.base))
  
  if(method == "auto"){
    # non-glycine whiskers ----
    if(x$glycine[1] == FALSE){
      # assign grouping for the regrowth section 
      # (identify section of the whisker corresponding to the deployment period since all whiskers have a pre-deployment section from the tip)
      # get row number corresponding to start.date by accounting for root length from the tip. 
      l <- x %>% 
        arrange(num.from.tip) %>% 
        mutate(len.from.tip = cumsum(sample.len)) %>% .$len.from.tip
      rn <- which.min(abs(x$root.length[1] * ratio - l)) + 1 # 0.8 is the ratio of pre-deployment root length to the post-deployment root length. 
      
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
  } else if(method == "rootlength"){
    # assign grouping for the regrowth section 
    # get row number corresponding to start.date by accounting for root length from the tip. 
    l <- x %>% 
      arrange(num.from.tip) %>% 
      mutate(len.from.tip = cumsum(sample.len)) %>% .$len.from.tip
    rn <- which.min(abs(x$root.length[1] * ratio - l)) + 1
    
    x1 <- x %>% 
      arrange(num.from.tip) %>% 
      mutate(rg = 0)
    x1$rg[rn:nrow(x1)] <- seq(1,nrow(x1)-rn+1,1)
    
    # calculate growth.rate and resolution ---
    if(x$plucked[1] == TRUE){
      # plucked whiskers
      x <- x1 %>% mutate(res = round(difftime(end.date,start.date, units = 'days')/(original.length - root.length * ratio), 2),
                         growth.rate = (original.length - root.length)/as.integer(difftime(end.date, start.date, units = 'days')))
    } else {
      # non-glycine whiskers
      x <- x1 %>% mutate(res = round(difftime(end.date,start.date, units = 'days')/original.length, 2),
                         growth.rate = original.length/as.integer(difftime(end.date, start.date, units = 'days')))
    }  
  }
    
  # estimate dates ------------------------------------------
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
  x2 <- full_join(x, x1,
                  by = c("num.from.base", "c13", "n15", "sample.len", "root", "root.length", "original.length", "plucked", "start.date", "end.date", "glycine", "has.glycine", "num.from.tip", "rg", "res", "growth.rate")) %>% 
    filter(is.na(date.est))
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
  mutate(data = purrr::map(data, interp.regrowth.dates, method = 'rootlength', ratio = 0.8)) %>% 
  unnest()

# look at growth rates 
growth <- d2 %>% 
  group_by(seal) %>% 
  summarise(growth.rate = first(growth.rate), glycine = first(glycine)) 
# growth
growth %>% group_by(glycine) %>% 
  summarise(mean = mean(growth.rate), sd = sd(growth.rate)) # for glycine seals, growth rate using glycine market > using root length method

# check if there is delay in glycine assimilation to appear on whisker. if method = 'root length" 
d2 %>% filter(has.glycine == TRUE) %>%
  group_by(seal) %>%
  summarise(injection.date = first(start.date),
            estimated.start.date = min(date.est),
            delta.dates = difftime(estimated.start.date, injection.date))




# Plot spatial map of stable isotope values -------------------------------
# Load time spent in cell data from home range analyses
load("~/Dropbox/phd/LNFS research/LNFS processed data/animal home range/time spent/timespent_and_percentage_by_seal_season_trip.RData")

# assign week and month to each trip
h <- d %>% 
  group_by(seal, trip) %>% 
  mutate(month = month(median(c(startdate, enddate)), label = T, abbr = T),
         week = week(median(c(startdate, enddate))),
         median.date = as.Date(median(c(startdate,enddate))), 
         q = quantile(timespent, .9)) %>% 
  filter(seal %in% unique(d2$seal))

# get nearest matching track date for SIA
d3 <- d2
d3$n15[d3$has.glycine == TRUE | d3$date.est < d3$start.date] <- NA # remove glycine values and pre-deployment values
d3$c13[ d3$date.est < d3$start.date] <- NA # remove glycine values and pre-deployment values
d3 <- d3 %>% group_by(seal) %>% filter(num.from.base != 1)
data <- split(d3, paste0(d3$seal, d3$num.from.base))
d3 <- purrr::map(data, function(x){
  id <- x$seal
  track <- h %>% ungroup() %>%  dplyr::filter(seal == id)
  dates <- unique(track$median.date)
  nearest_date <- dates[which.min(abs(x$date.est - dates))]
  if(length(nearest_date) > 0) x$median.date <- nearest_date
  x
}) %>% reduce(bind_rows)

d3 <- d3 %>% 
  filter(!is.na(median.date)) %>% 
  group_by(seal, median.date) %>% 
  summarise(c13 = mean(c13, na.rm = T), n15 = mean(n15, na.rm = T))


# home range to SIA by week
h1 <- left_join(h, d3)
h2 <- h1 %>% group_by(x, y, year, month) %>% 
  filter( timespent >= q) %>% 
  summarise(c13 = mean(c13, na.rm = T), n15 = mean(n15, na.rm = T))

# plot
pdf(file = './plots/SI_spatialmaps.pdf', width = 10, height = 7)
library(marmap)
map <- map_data("worldHires")
bathy <- marmap::getNOAA.bathy(129,145,-43, -35,res=1, keep=TRUE)
bathy_map <- tbl_df(fortify(bathy))
b2000 <- bathy_map %>% filter(z >= -2000 & z <= -1000) %>% group_by(x) %>% summarise(y = mean(y))

mytheme <- theme(legend.position = 'bottom',
                 text = element_text(size = 8.5),
                 panel.background = element_rect(fill = "white", colour = 'black'),
                 panel.border = element_rect(linetype = "solid", fill = NA))
p1 <- h2 %>% ungroup() %>% 
  filter(!is.na(n15)) %>% 
  ggplot(aes(x, y))  +
  geom_raster(aes(fill = n15)) + 
  geom_map(map_id = "Australia", map = map, fill = 'grey70') +
  geom_path(data = b2000, colour = 'grey70') +
  facet_wrap(year ~ month) + 
  scale_fill_viridis(option = 'C') + 
  mytheme

p2 <- h2 %>% ungroup() %>% 
  filter(!is.na(c13)) %>% 
  ggplot(aes(x, y))  +
  geom_raster(aes(fill = c13)) + 
  geom_map(map_id = "Australia", map = map, fill = 'grey70') +
  geom_path(data = b2000, colour = 'grey70') +
  facet_wrap(year ~ month) + 
  scale_fill_viridis(option = 'C') + 
  mytheme


p3 <- ggarrange(p1, p2, ncol = 2)
print(p3)

seals <- unique(d2$seal[d2$seal %in% h1$seal])
for(i in seals){
  p1 <- h1 %>% ungroup() %>% 
    filter(!is.na(n15), seal == i) %>% 
    filter( timespent >= q) %>% 
    ggplot(aes(x, y))  +
    geom_raster(aes(fill = n15)) + 
    geom_map(map_id = "Australia", map = map, fill = 'grey70') +
    geom_path(data = b2000, colour = 'grey70') +
    facet_wrap(seal ~ month) + 
    scale_fill_viridis(option = 'C') + 
    mytheme + 
    lims(x = c(range(h2$x)), y = c(range(h2$y)))
  
  p2 <- h1 %>% ungroup() %>% 
    filter(!is.na(c13), seal == i) %>% 
    filter( timespent >= q) %>% 
    ggplot(aes(x, y))  +
    geom_raster(aes(fill = c13)) + 
    geom_map(map_id = "Australia", map = map, fill = 'grey70') +
    geom_path(data = b2000, colour = 'grey70') +
    facet_wrap(seal ~ month) + 
    scale_fill_viridis(option = 'C') + 
    mytheme + 
    lims(x = c(range(h2$x)), y = c(range(h2$y)))
  
  p3 <- ggarrange(p1, p2, ncol = 2)
  print(p3)
  print(i)
}

dev.off()

h1 %>% ggplot(aes(y, c13, group = y)) + 
  geom_boxplot() + facet_wrap(~year)

# Plots -------------------------------------------------------------------


# create sequence from tip (oldest info) 
d3 <- d2 %>% 
  arrange(seal, num.from.base) %>% 
  group_by(seal) %>% 
  mutate(num.from.tip = rev(num.from.base))

# modify estimated dates for plotting purposes
d3 <- d3 %>% 
  mutate(date = as.Date(format(date.est, format = '2017-%m-%d'))) %>% 
  mutate(date = ifelse(date > ymd(20171101), 
                       format(date, format = '2016-%m-%d'), 
                       format(date, format = '2017-%m-%d'))) %>% 
  mutate(date = as.Date(date)) %>% 
  mutate(year = year(start.date)) %>% 
  mutate(start.date = as.Date(format(start.date, format = '2017-%m-%d')), end.date = as.Date(format(end.date, format = '2017-%m-%d')))
  
# import gls tracks for comparison
load("./data/glsTracksAll_cleaned.RData")
colnames(tracks) <- tolower(colnames(tracks))
tracks <- tracks %>% mutate(id = as.numeric(as.character(id)))

# plot c13 and n15 time series
# pdf(file = './plots/SI_plots.pdf', width = 10, height = 10)
pal <- wesanderson::wes_palette('BottleRocket2', 3, type = 'continuous' )
t1 <- min(d3$start.date) 
t2 <- max(d3$end.date)

# d13c vs d15n
d3 %>% 
  filter(has.glycine == FALSE) %>% 
  ggplot(aes(x = c13, y = n15, colour = num.from.tip)) +
  geom_point(aes(size = num.from.tip)) +
  facet_wrap(~seal) + 
  scale_colour_viridis(option = 'B') + 
  # scale_colour_gradient2(low = pal[1], mid = pal[2], high = pal[3], midpoint = 13) +
  theme_bw()

# date vs d13c
d3 %>% 
  ggplot(aes(x = date, y = c13)) + 
  geom_point(size = 1, aes(colour = has.glycine)) + 
  geom_line(size = 0.5) +
  geom_vline(aes(xintercept = start.date), linetype = 'dashed') +
  geom_vline(aes(xintercept = end.date), linetype = 'dashed') +
  labs(x = 'date') +
  scale_x_date(date_breaks = '1 month', date_labels = '%b', limits = c(t1, t2)) +
  theme_bw() + 
  facet_wrap(~seal)

# date vs d15n (without glycine segments)
d3 %>% 
  filter(has.glycine == FALSE) %>% 
  ggplot(aes(x = date, y = n15)) + 
  geom_point(size = 1) + 
  geom_line(size = 0.5) +
  geom_vline(aes(xintercept = start.date)) +
  geom_vline(aes(xintercept = end.date)) +
  labs(x = 'date') +
  scale_x_date(date_breaks = '1 month', date_labels = '%b', limits = c(t1,t2)) +
  theme_bw() + 
  facet_wrap(~seal)

# date vs d15n (only glycine seals)
d3 %>% 
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

pdf(file = './plots/individual_SI_plots_rootlengthMethod.pdf', width = 10, height = 7)
# each seal c13 and n15 with lat, lon
seals <- unique(d3$seal)[unique(d3$seal) %in% tracks$id]
for(i in seq_along(seals)){
  x <- seals[i]
  p1 <- d3 %>% filter(seal == x) %>% 
    ggplot(aes(x = date, y = c13)) + 
    geom_point(size = 1, aes(colour = has.glycine)) + 
    geom_line(size = 0.5) +
    geom_vline(aes(xintercept = start.date), linetype = 'dashed') +
    geom_vline(aes(xintercept = end.date), linetype = 'dashed') +
    labs(x = 'date', title = x) +
    scale_x_date(date_breaks = '1 month', date_labels = '%b', limits = c(t1, t2)) +
    theme_bw() + 
    theme(legend.position = c(0.8, 0.25),
          legend.background = element_rect(fill = alpha('white', 0), colour = NA),
          # legend.key.size = unit(.25,'cm'),
          legend.text = element_text(size = rel(1)),
          legend.title = element_text(size = rel(1)),
          legend.direction = 'horizontal',
          legend.spacing = unit(0, 'lines'),
          legend.margin = margin(3,3,3,3))
  
  p2 <- d3 %>% filter(seal == x) %>% 
    filter(has.glycine == FALSE) %>% 
    ggplot(aes(x = date, y = n15)) + 
    geom_point(size = 1, aes(colour = has.glycine)) + 
    geom_line(size = 0.5) +
    geom_vline(aes(xintercept = start.date), linetype = 'dashed') +
    geom_vline(aes(xintercept = end.date), linetype = 'dashed') +
    scale_x_date(date_breaks = '1 month', date_labels = '%b', limits = c(t1, t2)) +
    theme_bw() + 
    theme(legend.position = c(0.8, 0.25),
          legend.background = element_rect(fill = alpha('white', 0), colour = NA),
          # legend.key.size = unit(.25,'cm'),
          legend.text = element_text(size = rel(1)),
          legend.title = element_text(size = rel(1)),
          legend.direction = 'horizontal',
          legend.spacing = unit(0, 'lines'),
          legend.margin = margin(3,3,3,3),
          axis.title.x = element_blank())
  
  p3 <- tracks %>% filter(id == x) %>% 
    mutate(date = as.Date(format(date, format = '2017-%m-%d %H:%M:%S'))) %>% 
    ggplot(aes(x = date, y = lat)) + 
    geom_point(size = 0.3) + 
    geom_line(size = 0.5) +
    labs(x = 'date') +
    scale_x_date(date_breaks = '1 month', date_labels = '%b', limits = c(t1, t2)) +
    theme_bw() + 
    theme(axis.title.x = element_blank())
  
  p4 <- tracks %>% filter(id == x) %>% 
    mutate(date = as.Date(format(date, format = '2017-%m-%d %H:%M:%S'))) %>% 
    ggplot(aes(x = date, y = lon)) + 
    geom_point(size = 0.3) + 
    geom_line(size = 0.5) +
    scale_x_date(date_breaks = '1 month', date_labels = '%b', limits = c(t1, t2)) +
    theme_bw() + 
    theme(axis.title.x = element_blank())
  
  
  p5 <- tracks %>% filter(id == x) %>% 
    mutate(month = month(date)) %>% 
    ggplot(aes(lon, lat)) + 
    geom_path() + 
    facet_wrap(~month) + 
    geom_map(map_id = "Australia", map = map, fill = 'grey30') + 
    theme_bw() + 
    lims(x = c(range(tracks$lon)), y = c(range(tracks$lat)))
    
  
  tmp <- ggarrange(p1,p2,p3,p4, align = 'hv', nrow = 4)
  out <- ggarrange(tmp, p5, ncol = 2, widths = c(1,2)) 
  print(out)
  print(i)
}

dev.off()

# make a monthly tracks plot with all individuals
tracks %>%
  mutate(month = month(date)) %>% 
  ggplot(aes(lon, lat)) + 
  geom_path(aes(colour =  as.factor(id))) + 
  facet_wrap(~month) + 
  geom_map(map_id = "Australia", map = map, fill = 'grey30') + 
  theme_bw() + 
  lims(x = c(range(tracks$lon)), y = c(range(tracks$lat)))

# outliers -------------------------------------------------------------------
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

# pups sia ----------------------------------------------------------------
d0 <- read_csv("./data/lnfs_all_whisker_bulk_sia_correctedsequence.csv")
p <- d0 %>% filter(adult == FALSE) %>% 
  select(seal, num.from.base, c13, n15, sample.len, root, root.length, original.length, plucked, start.date, end.date, glycine, growth.rate, res) %>% 
  mutate(start.date = as.Date(strptime(start.date, format = '%d/%m/%y')), end.date = as.Date(strptime(end.date, format = '%d/%m/%y'))) %>% 
  arrange(seal, num.from.base)

m <- d2 %>% filter(seal %in% c(71,73)) # pup's mothers

p <- p %>% group_by(seal) %>% 
    mutate(cum.len = cumsum(sample.len)) %>% 
  mutate(date.est = end.date - res*lag(cum.len))
p$date.est[is.na(p$date.est)] <-p$end.date[is.na(p$date.est)] 

p %>% 
  filter(c13 > -18) %>%
  ggplot(aes(x = date.est, y = c13, colour = seal)) + 
  geom_point() + 
  geom_line() + 
  scale_x_date(date_breaks = '1 month', date_labels = '%b')
