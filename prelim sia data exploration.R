## explore whisker sia data ##
rm(list = ls())
library(tidyverse)
library(lubridate)
# library(gridExtra)

# load data
# note: for preliminary seals, pups were sectioned from root to tip, whereas adults were sectioned from tip (oldest) to root (newest).
# but sequence is already corrected mannually in csv from root (num 1) to tip
# adults analysed were regrowths
# pup whiskers were plucked at the end of fieldwork
d0 <- read_csv("./preliminary animals/lnfs_prelim_whisker_bulk_sia_correctedsequence.csv")
d1 <- d0 %>% 
  select(seal, num, c13, n15, sample.len, date) %>% 
  mutate(date = as.Date(strptime(date, format = '%d/%m/%y'))) %>% 
  arrange(seal, num)

# add rows for missing segments
# d1 <- lapply(split(d, d$seal), function(x) {
#   merge(expand.grid(seal = unique(x$seal), num = min(x$num):max(x$num)), x, all=TRUE)
#   })

# detect outliers
d1 %>% ggplot(aes(x = seal, y = c13)) +
  geom_boxplot() # pup 103 c13 obvious outlier 
d1 %>% ggplot(aes(x = seal, y = n15)) +
  geom_boxplot()

x <- d1$c13[d1$seal == '103']
boxplot.stats(x)$out
outlier_values <- boxplot.stats(x)$out  # outlier values.
boxplot(x, boxwex=0.1)
mtext(paste("Outliers: ", paste(outlier_values, collapse=", ")), cex=0.6)

# remove outliers
d1$c13[d1$seal == '103'][which(d1$c13[d1$seal == '103'] %in% outlier_values )] <- NA

# create groups for mum-pup pairs
d1$pair <- 1
d1$pair[d1$seal == '104' | d1$seal == '77'] <- 2
d1$pair[d1$seal == '103' | d1$seal == '71'] <- 3
d1$pair <- as.factor(d1$pair)

# interpolate dates of whisker segments ---
# reverse sequence from oldest (tip) to newest (root)
d1 <- d1 %>% 
  arrange(seal, desc(num)) %>% 
  group_by(seal) %>% 
  mutate(cum.sample.len = cumsum(sample.len))

# for adult regrowths
interp.regrowth.dates <- function(df){
  x <- df
  start.date <- min(x$date, na.rm = T)
  end.date <- max(x$date, na.rm = T)
  growth.rate <- difftime(end.date, start.date, units = 'days')/max(x$cum.sample.len)
  interp.dates <- start.date + growth.rate * lag(x$cum.sample.len)
  interp.dates <- c(start.date, na.omit(interp.dates))
  x$date.interp <- interp.dates
  x$growth.rate <- growth.rate
  x
}
d2 <- d1 %>% 
  group_by(seal) %>% 
  filter(!seal %in% c('103', '104')) %>% 
  nest() %>% 
  mutate(data = map(data, interp.regrowth.dates)) %>% 
  unnest()

d2 %>% group_by(seal) %>% summarise(growth.rate = first(growth.rate))

# for pup plucked whiskers
d3 <- d1 %>% 
  filter(seal %in% c('103', '104')) %>% 
  mutate(growth.rate = ifelse(seal == '103', 2.62, 3.07)) # manually input growth rate (calculate from regrowths)

interp.whisker.dates <- function(x){
  end.date <- max(x$date, na.rm = T)
  growth.rate <- x$growth.rate[1]
  start.date <- end.date - max(x$cum.sample.len) * growth.rate
  interp.dates <- start.date + growth.rate * lag(x$cum.sample.len)
  interp.dates <- c(start.date, na.omit(interp.dates))
  x$date.interp <- interp.dates
  x
}

d3 <- d3 %>% 
  group_by(seal) %>% 
  nest() %>% 
  mutate(data = map(data, interp.whisker.dates)) %>% 
  unnest()

# combine adults and pups
d1 <- bind_rows(d2,d3)


# plot c13 and n15 time series
d1 %>% 
  ggplot(aes(x = date.interp, y = c13)) + 
  geom_point(aes(color = seal), size = 1) + 
  geom_line(aes(color = seal), size = 0.5, alpha = 0.8) +
  labs(x = 'date') +
  scale_x_date(date_breaks = '1 month', date_labels = '%b') +
  theme_bw() +
  facet_wrap(~pair)

d1 %>% 
  ggplot(aes(x = date.interp, y = n15)) + 
  geom_point(aes(color = seal), size = 1) + 
  geom_line(aes(color = seal), size = 0.5, alpha = 0.8) +
  labs(x = 'date') +
  scale_x_date(date_breaks = '1 month', date_labels = '%b') +
  theme_bw() +
  facet_wrap(~pair)

d1 %>% 
  filter(seal %in% c('71', '77', '450', '450b')) %>% 
  ggplot(aes(x = date.interp, y = c13, colour = seal)) + 
  geom_point(size = 0.5) + 
  geom_line() +
  scale_x_date(date_breaks = '1 month', date_labels = '%b') +
  theme_bw()

pdf(file = 'lnfs_preliminary_SIA.pdf', width = 7, height = 3.5)
d1 %>% filter(!seal  %in% c('103', '104')) %>% 
  ggplot(aes(num, n15)) + 
  geom_line() + 
  geom_point(size = 1) +
  facet_wrap(~seal)

d1 %>% filter(!seal  %in% c('103', '104')) %>% 
  ggplot(aes(num, c13)) + 
  geom_line() + 
  geom_point(size = 1) +
  facet_wrap(~seal)
dev.off()


# plot
p1 <- ggplot(d, aes(x = cum_len, y = d13C)) + 
  geom_point(aes(color = pair, shape = seal), size = 1.3) + 
  geom_line(aes(color = pair, group = seal, linetype = seal), size = 0.5, alpha = 0.8) +
  labs(x = "whisker length tip to base (mm)")
  
p2 <- ggplot(d, aes(x = cum_len, y = d15N)) + 
  geom_point(aes(color = pair, shape = seal), size = 1.3) + 
  geom_line(aes(color = pair, group = seal, linetype = seal), size = 0.5, alpha = 0.8) +
  labs(x = "whisker length tip to base (mm)")

p3 <- ggplot(d, aes(x = d13C, y = d15N)) + 
  geom_point(aes(color = pair, shape = seal), size = 2)

# p3 <- ggplot(d, aes(x = cum_len, y = diffC)) + 
#   geom_point(aes(color = pair, shape = seal), size = 0.5) + 
#   geom_line(aes(color = pair)) + 
#   labs(x = "whisker length", y = 'd13C diff between current and previous segment')
# 
# p4 <- ggplot(d, aes(x = cum_len, y = diffN)) + 
#   geom_point(aes(color = pair), size = 0.5) + 
#   geom_line(aes(color = pair, linetype = seal)) + 
#   labs(x = "whisker length", y = 'd15N diff between current and previous segment')

tiff(filename = 'preliminary whisker SIA d13C d15N mum-pup pairs.tiff', width = 7.5, height = 7.5, units="in", res=300)
grid.arrange(p1, p2, ncol = 1, nrow = 2)
dev.off()

tiff(filename = 'preliminary whisker SIA d13C x d15N.tiff', width = 7.5, height = 7.5, units="in", res=300)
p3 
dev.off()
