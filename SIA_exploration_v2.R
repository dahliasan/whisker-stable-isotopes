## explore whisker sia data ##

# Load libraries ----------------------------------------------------------


rm(list = ls())
library(tidyverse)
library(lubridate)
library(ggpubr)
library(viridisLite)
library(mapdata)
library(viridis)
library(devtools)
# devtools::install_github('cttobin/ggthemr')
library(ggthemr)
options(dplyr.width = Inf) #enables head() to display all coloums

## Load data
# adults analysed were regrowths
# pup whiskers were plucked at the end of fieldwork
d0 <- read_csv("./data/lnfs_all_whisker_bulk_sia_correctedsequence.csv")
d1 <- d0 %>% 
  filter(adult == TRUE) %>% 
  select(seal, num.from.base, c13, n15, sample.len, root, root.length, original.length, plucked, start.date, end.date, glycine) %>% 
  mutate(start.date = as.Date(strptime(start.date, format = '%d/%m/%y')), end.date = as.Date(strptime(end.date, format = '%d/%m/%y'))) %>% 
  mutate(root = ifelse(is.na(root), FALSE, TRUE)) %>% 
  arrange(seal, num.from.base) %>% 
  mutate(year = year(start.date))

# Fix up data errors ------------------------------------------------------
## Remove seal 450, keep seal , rename 450b to 450
d1 <- d1 %>% filter(!seal == '450') 
d1$seal[d1$seal == '450b'] <- '450'

## Check for data errors 
d1 %>% filter(num.from.base > 0) %>% 
  group_by(seal) %>% 
  summarise(tsamplelen = sum(sample.len), originallen = first(original.length), diff = tsamplelen - originallen)

d1 %>% filter(plucked == T, root == T) %>% 
  group_by(seal) %>% 
  summarise(tsamplelen = sum(sample.len), originallen = first(root.length), diff = tsamplelen - originallen)


## Recalculate each sample length 
# 2017 seals ignore root section
d11 <- d1 %>% filter(root == FALSE, year == 2017) %>% 
  mutate(norootlen = original.length - root.length) %>% 
  group_by(seal) %>% 
  mutate(sample.len = norootlen/n())

# 2016 seals ignore root section
d12 <- d1 %>% filter(root == FALSE, year == 2016) %>% 
  mutate(norootlen = original.length) %>% 
  group_by(seal) %>% 
  mutate(sample.len = norootlen/n())

# 2016 and 2017 root sections
d13 <-  d1 %>% filter(root == T) %>% 
   group_by(seal) %>% 
   mutate(sample.len = root.length/n())

# join everything back into single dataframe - data errors fixed. 
d1 <- full_join(d11, d12) %>% full_join(d13) %>% 
  arrange(seal, num.from.base) %>% 
  select(-norootlen)

# Glycine cutoff ----------------------------------------------------------
## Find cutoff n15 value for glycine seals
d1 %>% 
  filter(glycine == TRUE) %>% 
  ggplot(aes(y = n15, x = num.from.base)) + geom_point(size = 0.5) + geom_line() + facet_wrap(~seal) + 
  geom_hline(aes(yintercept = 17.6))

## Label glycine segments
d1 <- d1 %>% mutate(has.glycine = ifelse(n15 > 17.6, TRUE, FALSE))
d1$has.glycine[is.na(d1$has.glycine)] <- FALSE



# Remove Outliers ---------------------------------------------------------
# Remove suspicious outlier data points in #319 
d1$c13[d1$seal == '319' & d1$n15 > 17 & d1$has.glycine == FALSE] <- NA
d1$n15[d1$seal == '319' & d1$n15 > 17 & d1$has.glycine == FALSE] <- NA


# Data exploration --------------------------------------------------------
d1 %>% filter(n15 <= 17.6) %>% 
  ggplot(aes(y = n15, x = seal, group = seal, colour = glycine)) + 
  geom_boxplot()

d1 %>% 
  ggplot(aes(y = c13, x = seal, group = seal, colour = glycine)) + 
  geom_boxplot()

d1 %>% 
  group_by(glycine) %>% 
  summarise(maxN = max(n15, na.rm = T)) # max d15N is 16 for KI shelf area (agrees with Lowther et al. 2013)

d1 %>% filter(has.glycine == FALSE) %>% 
  ggplot(aes(x = c13, y = n15, colour = glycine)) + 
  geom_point(size = 0.5) +
  facet_wrap(~seal)


## Between year differences
d1 %>% 
  filter(!seal %in% c('307', '317', '340')) %>% 
  ggplot(aes(year, c13, group = year)) + 
  geom_boxplot() + 
  labs(title = 'Interannual differences of lactating seals')

d1 %>% 
  filter(has.glycine == FALSE) %>% 
  filter(!seal %in% c('307', '317', '340')) %>% 
  ggplot(aes(year, n15, group = year)) + 
  geom_boxplot() + 
  labs(title = 'Interannual differences of lactating seals')
  

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
g1 <- d1 %>%
  filter(glycine == TRUE) 

# plot n15
g1 %>%
  ggplot(aes(num.from.base, n15)) + 
  geom_point(aes(colour = has.glycine)) + 
  geom_line() + 
  facet_wrap(~seal)

## Is the predeployment whisker lengthequal to the root length? ##
# i.e. length from tip to glycine start date = root length?
tmp <- g1 %>% 
  group_by(seal) %>% 
  filter(has.glycine == TRUE) %>% 
  mutate(gly.start = max(num.from.base)) %>% 
  select(seal, gly.start)
tmp <- tmp[!duplicated(tmp),]

g1 <- left_join(g1, tmp)

a <- g1 %>% 
  filter(num.from.base > gly.start) %>%
  filter(!seal %in% c('318', '322', '324')) %>%  # remove seals with missing predeployment sections
  group_by(seal) %>% 
  summarise(predeploy.len = sum(sample.len), 
            root.length = first(root.length), 
            diff = predeploy.len - root.length, 
            ratio = predeploy.len/root.length)

a %>% filter(seal != '326') %>% 
  summarise(mean.ratio = mean(ratio)) # pre-deployment root length average 0.8 of current root length. 


# Estimate dates of whisker segments --------------------------------------
# for adult regrowths
source('./estRegrowthDates_v1.R')

d2 <- d1 %>% 
  group_by(seal) %>% 
  nest() %>% 
  mutate(data = purrr::map(data, estRegrowthDates, method = 'rootlength', ratio = 1, estimateBackwards = TRUE)) %>% 
  unnest()

# look at growth rates 
growth <- d2 %>% 
  group_by(seal) %>% 
  summarise(growth.rate = first(growth.rate), glycine = first(glycine)) 
growth
growth %>% group_by(glycine) %>% 
  summarise(mean = mean(growth.rate), sd = sd(growth.rate)) # for glycine seals, growth rate using glycine marker > using root length method
min(growth$growth.rate)
# check if there is delay in glycine assimilation to appear on whisker. if method = 'rootlength" 

a <- d2 %>% 
  filter(!seal %in% c('322', '324')) %>% # these are seals with incomplete glycine segments. 
  filter(has.glycine == TRUE) %>% 
  group_by(seal) %>%
  summarise(injection.date = first(start.date),
            estimated.start.date = min(date.est),
            delta.dates = difftime(estimated.start.date, injection.date))
a
mean(a$delta.dates)
sd(a$delta.dates)


# Match SIA to Core foraging areas -------------------------------
# Load time spent in cell data from home range analyses
load("~/Dropbox/phd/LNFS research/LNFS processed data/animal home range/time spent/timespent_and_percentage_by_seal_season_trip.RData")

h <- d %>% 
  group_by(seal, trip) %>% 
  # Assign the month and calculate the 90th quantile of time spent for each trip
  mutate(median.date = as.Date(median(c(startdate,enddate))), 
         month = month(median.date, label = T, abbr = T),
         q = quantile(timespent, .9)) %>% 
  # subset seals with SIA data
  filter(seal %in% unique(d2$seal))

## Match SIA values to nearest foraging trip
# fix seal names
i <- which(str_count(d2$seal) == 2)
d2$seal[i] <- paste0('0', d2$seal[i])

d3 <- d2 %>% 
  group_by(seal) %>% 
  # remove segments closest to root (potentially anomalous values due to root effect)
  filter(num.from.base != min(num.from.base)) %>% 
  ungroup() 

# remove glycine values and pre-deployment values
d3$n15[d3$has.glycine == TRUE | d3$date.est < d3$start.date] <- NA 
d3$c13[ d3$date.est < d3$start.date] <- NA 

# for each whisker segment, find nearest foraging trip
data <- split(d3, paste0(d3$seal, d3$num.from.base))
d3 <- purrr::map(data, function(x){
  id <- x$seal
  track <- h %>% ungroup() %>%  dplyr::filter(seal == id)
  dates <- unique(track$median.date)
  nearest_date <- dates[which.min(abs(x$date.est - dates))]
  if(length(nearest_date) > 0) x$median.date <- nearest_date
  return(x)
}) %>% reduce(bind_rows)


# Plot study region map ---------------------------------------------------
## Plot tracks showing core foraging areas
library(marmap)
map <- map_data("worldHires", "Australia")
bathy <- marmap::getNOAA.bathy(129,151,-46, -35,res=1, keep=TRUE)
bathy_map <- tbl_df(fortify(bathy))
b2000 <- bathy_map %>% filter(z >= -2000 & z <= -1000) %>% group_by(x) %>% summarise(y = mean(y))

mytheme <- theme(legend.position = 'top',
                 text = element_text(size = 10),
                 panel.background = element_rect(fill = "white", colour = 'black'),
                 panel.border = element_rect(linetype = "solid", fill = NA))

ggthemr("fresh")
a <- h %>% 
  ungroup() %>% 
  mutate(trip = as.numeric(trip)) %>% 
  group_by(seal, trip) %>% 
  mutate(core = ifelse(timespent >= q & y <= min(y) + 2, 'core', 'non-core'),
         seal_trip = paste(seal, trip, sep = '_')) %>% 
  filter(seal == 450) %>% 
  ggplot(aes(x, y)) + 
  geom_path(data = b2000, colour = 'grey70') +
  geom_raster(aes(fill = core)) + 
  geom_map(map_id = "Australia", map = map, fill = 'grey70') +
  facet_wrap(~trip) + 
  mytheme + 
  lims(x = c(134, 143), y = c(-43, - 36)) + 
  labs(x = 'Lon', y = "Lat", fill = 'Foraging area type') 
# ggthemr_reset()

b <- ggplot() +
  geom_map(data = map, aes(x = long, y = lat, map_id = region),
           map = map, colour = NA, fill = "grey60") +
  geom_rect(data = data.frame(),
            aes(xmin = 134, xmax = 143, ymin = -43, ymax = -36),
            colour = "red", fill = NA) +
  geom_text(aes(label = 'GAB', x = 130, y = -33.2), colour = 'black', size = 2.5) + 
  coord_map(xlim = c(110, 155), ylim = c(-45, -10)) +
  labs(x = NULL, y = NULL) + 
  theme_pubclean()

maptheme <- theme(
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank(),
  panel.border = element_rect(fill = NA, colour = "black"),
  panel.background = element_blank()
)

library(grid)
grid.newpage()
vp_a <- viewport(width = 1, height = 1, x = 0.5, y = 0.5)  # the larger map
vp_b <- viewport(width = 0.2, height = 0.2, x = 0.8, y = 0.12)  # the inset in upper right
print(a, vp = vp_a)
print(b +maptheme , vp = vp_b)

# Plot spatial map --------------------------------------------------------
# Get only core foraging locations and summarise SI values by year and month
h1 <- h %>% 
  group_by(seal, trip, median.date) %>%
  filter(timespent >= q) %>% 
  filter(y <= min(y) + 2) %>%
  summarise_at(vars(x,y), funs(median(., na.rm = T)))
h2 <- left_join(d3, h1)
h2 <- mutate(h2, year = as.factor(year(start.date)), month = month(median.date, abbr = T, label = T))
h2 <- h2 %>% group_by(x, y, year, month) %>% 
  summarise(c13_sd = sd(c13, na.rm = T), n15_sd = sd(n15, na.rm = T),
            c13 = mean(c13, na.rm = T), n15 = mean(n15, na.rm = T))

# plot
# pdf(file = './plots/SI_spatialmaps_2.pdf', width = 10, height = 7)


# Overall mean spatial map for each year
p1 <- h2 %>% 
  ungroup() %>% 
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


p1 <- h2 %>% ungroup() %>% 
  filter(!is.na(n15_sd)) %>% 
  ggplot(aes(x, y))  +
  geom_raster(aes(fill = n15_sd)) + 
  geom_map(map_id = "Australia", map = map, fill = 'grey70') +
  geom_path(data = b2000, colour = 'grey70') +
  facet_wrap(year ~ month) + 
  scale_fill_viridis(option = 'C') + 
  mytheme


p2 <- h2 %>% ungroup() %>% 
  filter(!is.na(c13_sd)) %>% 
  ggplot(aes(x, y))  +
  geom_raster(aes(fill = c13_sd)) + 
  geom_map(map_id = "Australia", map = map, fill = 'grey70') +
  geom_path(data = b2000, colour = 'grey70') +
  facet_wrap(year ~ month) + 
  scale_fill_viridis(option = 'C') + 
  mytheme


p3 <- ggarrange(p1, p2, ncol = 2)
print(p3)


# Individual seal spatial map
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



# Latitudal pattern -------------------------------------------------------

h1 <- h %>% 
  group_by(seal, trip, median.date) %>%
  filter(timespent >= q) %>% 
  filter(y <= min(y) + 2) %>%
  summarise_at(vars(x,y), funs(median(., na.rm = T)))
h2 <- left_join(d3, h1)
h2$year <- as.factor(year(h2$start.date))

mytheme <- theme(legend.position = 'bottom',
                 text = element_text(size = 8.5),
                 panel.background = element_rect(fill = "white", colour = 'black'),
                 panel.border = element_rect(linetype = "solid", fill = NA))
# h2 %>% ungroup() %>% 
#   mutate(y2 = signif(y, 2)) %>%
#   filter(!is.na(n15)) %>% 
#   ggplot(aes(y = n15))  +
#   geom_boxplot(aes(x = abs(y2), group = y2)) + 
#   geom_smooth(aes(x = abs(y)), linetype = 'dashed', size = 0.5) + 
#   facet_wrap(~year) + 
#   mytheme + 
#     labs(x = 'latitude south')

h2 %>% ungroup() %>%
  filter(!is.na(x)) %>%
  mutate(y2 = signif(y, 2)) %>%
  ggplot(aes(y = n15))  +
  geom_boxplot(aes(x = abs(y2), group = y2)) +
  geom_smooth(aes(x = abs(y)), linetype = 'dashed', size = 0.5) +
  geom_point(aes(x = abs(y))) + 
  facet_wrap(~year) +
  mytheme +
    labs(x = 'latitude south')

h2 %>% ungroup() %>%
  filter(!is.na(x)) %>%
  mutate(y2 = signif(y, 2)) %>%
  ggplot(aes(y = c13))  +
  geom_boxplot(aes(x = abs(y2), group = y2)) +
  geom_smooth(aes(x = abs(y)), linetype = 'dashed', size = 0.5, method = 'lm') +
  geom_point(aes(x = abs(y))) + 
  facet_wrap(~year) +
  mytheme +
  labs(x = 'latitude south')


## Fit GAM 
library(mgcv)
lmc <- lmeControl(niterEM=300)
m1 <- gamm(n15 ~ s(y, by = year), random = list(seal = ~1), data = h2)
summary(m1$gam)
par(mfrow=c(2,2))
gam.check(m1$gam)
par(mfrow=c(1,1))
m2 <- gamm(abs(c13) ~ s(y, by = year), random = list(seal = ~1), data = h2,  family = Gamma(link = 'log'), control = lmc)
summary(m2$gam)
par(mfrow=c(2,2))
gam.check(m2$gam)
par(mfrow=c(1,1))

m3 <- gamm(abs(c13) ~ s(y, by = year), random = list(seal = ~1), data = h2,  family = Gamma(link = 'log'), control = lmc)
AIC(m2$lme,m3$lme)

# make plots
par(mfrow=c(2,2), mar=c(5,5,1,1))
plot(m2$gam, select = 1, xlab = NA, ylab = 's(lat,1):year2016' )
mtext(expression(paste("(a) ", delta ^"13", "C")), adj = 0, line = -2, at = c(-45, 0.5))
plot(m2$gam, select = 2, xlab = NA, ylab = 's(lat,1):year2017' )
mtext(expression(paste("(b) ", delta ^"13", "C")), adj = 0, line = -2, at = c(-45, 0.5))
plot(m1$gam, select = 2, xlab = 'Latitude (degrees)', ylab = 's(lat,3.58):year2017')
mtext(expression(paste("(c) ", delta ^"15", "N")), adj = 0, line = -2, at = c(-45, 0.5))


# Compare shelf vs oceanic SI ---------------------------------------------
# Extract TOPO for locations
library(raster)
library(marmap)
bathy <- getNOAA.bathy(129,151,-46, -35, res=1, keep=TRUE)
bathy_map <- tbl_df(fortify(bathy))
bathy_ras <- bathy_map
coordinates(bathy_ras) <- ~x + y
gridded(bathy_ras) <- TRUE
bathy_ras <- raster(bathy_ras)

h1 <- h %>% 
  group_by(seal, trip, median.date) %>%
  filter(timespent >= q) %>% 
  filter(y <= min(y) + 2) %>%
  summarise(x = median(x, na.rm = T), y = median(y, na.rm = T), t.timespent = sum(timespent))

h2 <- left_join(d3, h1)
h2$year <- as.factor(year(h2$start.date))
h2$z <- raster::extract(bathy_ras, data.frame(x = h2$x, y = h2$y))
# h2 %>% ggplot(aes(x, y)) + geom_raster(aes(fill = z)) + scale_fill_viridis()
h2 <- h2 %>% mutate(region = ifelse(z < -2000, 'oceanic', 'shelf'))

## Manually assign shelf/oceanic group to seal 318 and 329 - based on TDR data
h2$region[h2$seal == '318' & h2$date.est < ymd(20170601)] <- 'shelf'
h2$region[h2$seal == '329' & h2$date.est < ymd(20170601)] <- 'oceanic'

## Check if region assigned correctly: 
## Should mostly have 069 and 071 seals without region because they don't have tracks. 
h2 %>% filter(is.na(region)) %>% dplyr::select(seal, region) %>% View()

h2 %>% filter(!is.na(region)) %>% 
  ggplot(aes(x = region, y = n15,  group = region)) + geom_boxplot()
h2 %>% filter(!is.na(region)) %>% 
  ggplot(aes(x = region, y = c13,  group = region)) + geom_boxplot()

## T-test
# Check for normality
shapiro.test(h2$n15)
shapiro.test(h2$c13)

# Non-normal t-test
a <- h2 %>% filter(region == 'shelf') %>% pull(c13) %>% na.omit()
b <- h2 %>% filter(region == 'oceanic') %>% pull(c13) %>% na.omit()
wilcox.test(a,b) # p <0.05
c(mean(a), sd(a), mean(b), sd(b)) %>% round(2)


a <- h2 %>% filter(region == 'shelf') %>% pull(n15) %>% na.omit()
b <- h2 %>% filter(region == 'oceanic') %>% pull(n15) %>% na.omit()
wilcox.test(a,b) # p <0.05
c(mean(a), sd(a), mean(b), sd(b))



# SI vs Habitat + Year ----------------------------------------------------
# Instead of doing a buncha t-tests as before, use GLM instead
library(lme4)
library(effects)
library(bbmle)

## Including root sections in 2017 
glm_data <- h2 %>% 
  filter(!is.na(region)) %>% 
  mutate(region = as.factor(region))

## 13C 
m0 <- glm(c13 ~ region*year, data = glm_data)
m1<- lmer(c13 ~ region*year + (1|seal), data = glm_data)
anova(m1, m0) # yes keep seal random effect
drop1(m1, test="Chisq")
m2 <- update(m1, .~. - region:year)
drop1(m2, test = 'Chisq')
ICtab(m0,m1,m2, base = T, logLik = T, type = 'AIC') 
summary(m2)
plot(allEffects(m2))

ICtab(m0,m1,m2, base = T, logLik = T, type = 'AIC') 
m2 %>% fit2df()

## 15N
mm0 <- glm(n15 ~ region*year, data = glm_data)
mm1<- lmer(n15 ~ region*year + (1|seal), data = glm_data)
anova(mm1, mm0) # yes keep seal random effect
drop1(mm1, test="Chisq")
mm2 <- update(mm1, .~. - region:year)
ICtab(mm0,mm1,mm2, base = T, logLik = T, type = 'AIC') 
summary(mm1)
plot(allEffects(mm1))

write.csv(rbind(tidy(m2), tidy(mm1)), './output/GLM_final_fit.csv')

## Excluding root sections in 2017
glm_data <- h2 %>% 
  filter(!is.na(region), root == FALSE) %>% 
  mutate(region = as.factor(region))

## 13C 
m0 <- glm(c13 ~ region*year, data = glm_data)
m1<- lmer(c13 ~ region*year + (1|seal), data = glm_data)
anova(m1, m0) # yes keep seal random effect
drop1(m1, test="Chisq")
m2 <- update(m1, .~. - region:year)
drop1(m2, test = 'Chisq')
ICtab(m0,m1,m2, base = T, logLik = T, type = 'AIC') 
summary(m2)
plot(allEffects(m2))

ICtab(m0,m1,m2, base = T, logLik = T, type = 'AIC') 
m2 %>% fit2df()

## 15N
mm0 <- glm(n15 ~ region*year, data = glm_data)
mm1<- lmer(n15 ~ region*year + (1|seal), data = glm_data)
anova(mm1, mm0) # yes keep seal random effect
drop1(mm1, test="Chisq")
mm2 <- update(mm1, .~. - region:year)
ICtab(mm0,mm1,mm2, base = T, logLik = T, type = 'AIC') 
summary(mm1)
plot(allEffects(mm1))

ICtab(mm0,mm1,mm2, base = T, logLik = T, type = 'AIC') 

write.csv(rbind(tidy(m2), tidy(mm1)), './output/GLM_final_fit_exclude_root.csv')


# Summarise for manuscript ------------------------------------------------
summ <- d2 %>% group_by(seal) %>% 
  summarise(startdate = first(start.date),
            enddate = first(end.date),
            glycine = first(glycine), 
            growthrate = first(growth.rate),
            tlength = first(original.length),
            rlength = first(root.length),
            nsections = n())

a <- summ %>% pull(growthrate) # non-glycine method
b <- c(0.186, 0.271, 0.136, 0.246, 0.186, 0.151, 0.155, 0.17)

wilcox.test(a,b) # p  > 0.05
round(c(mean(a), sd(a), mean(b), sd(b)),2)

# MCLUST ------------------------------------------------------------------
library(mclust)
library(ellipse)
# Fit the model
data <- d3
data <- data %>% filter(!is.na(n15), !is.na(c13))
seals <- data %>% pull(seal)
num.from.tip <- data %>% pull(num.from.tip)
data <- data %>% dplyr::select(c13,n15)
BIC <- mclustBIC(data)
plot(BIC)
summary(BIC)
m1 <- Mclust(data, x = BIC)
plot(m1, what = 'classification') 
data$classification = as.factor(m1$classification)
data$seal <- as.factor(seals)
data$num.from.tip <- num.from.tip

get.ellipses <- function(coords, mclust.fit){
  centers <- mclust.fit$parameters$mean[coords, ]
  vars <- mclust.fit$parameters$variance$sigma[coords, coords, ]
  purrr::map(1:ncol(centers), function(cluster){
    data.frame(ellipse(vars[,,cluster], centre = centers[, cluster],
                       level = 0.5), classification = cluster)
  })
}

mclust.el <- get.ellipses(c("c13", "n15"), m1) %>% 
  reduce(bind_rows) %>% 
  mutate(classification = as.factor(classification))

## Run compare shelf vs oceanic SI first to get h2
data <- left_join(data, h2) 
data$fs <- 'shelf/oceanic'
data$fs[data$seal %in% c('329', '353')] <- 'oceanic'
data$fs[data$seal %in% c('318', '315')] <- 'shelf'
data$fs[data$seal %in% c('307', '324', '317', '340')] <- 'unconstrained'
data <- data %>% mutate(year = as.factor(year(start.date)))
data$uncertainty <- m1$uncertainty

data$region[is.na(data$region)] <- "NA"
data %>%
  filter(!is.na(region)) %>%
  ggplot(aes(c13, n15)) +
  # geom_point(aes(shape = classification, colour = year), size = 2)+
  geom_point(aes(shape = classification, colour = region, size = uncertainty))+
  # geom_point(aes(shape = classification, colour = fs), size = 2)+
  geom_path(data = mclust.el, aes(group = classification, linetype = classification)) + 
  theme_pubr(base_size = 13, legend = "right") +  
  # scale_color_brewer(type = 'qual', palette = 7, direction = -1) + 
  scale_colour_manual(values = c('#ffa600', '#114d78', '#ce5291')) + 
  labs(x = expression(paste(delta ^"13", "C")), y = expression(paste(delta ^"15", "N")))


  
## Uncertainty
q <- quantile(m1$uncertainty)
uncertainty <- tbl_df(data.frame(uncertainty = m1$uncertainty, quantile = NA))
uncertainty$quantile <- cut(uncertainty$uncertainty, q, labels = c('q25', 'q50', 'q75', 'q100'))
uncertainty$quantile[which(is.na(uncertainty$quantile))] <- 'q25'
m1$uncertainty %>% hist()
plot(m1, what = 'uncertainty')

## Investigate different groups in oceanic habitat
## oceanic habitat classification are 1 and 3. 

# Load female size data 
size <- size %>% dplyr::select(animal_id_1, deploy_mass_kg, deploy_length_cm, deploy_girth_cm) %>% 
  rename(seal = animal_id_1, mass = deploy_mass_kg, length = deploy_length_cm, girth = deploy_girth_cm)

# Fix seal names
i <- which(str_count(size$seal) == 2)
size$seal[i] <- paste0('0', size$seal[i])

# add size data to model data
data_ocean <- data %>% 
  filter(!is.na(x)) %>%
  left_join(size) %>% 
  dplyr::select(seal, classification, year, length, x, y)

# GLMER doesn't work well... model nearly unidentifiable

## T-test

# Non-normal t-test
a <- data_ocean %>% filter(classification == '1') %>% pull(y) %>% na.omit()
b <- data_ocean %>% filter(classification == '3') %>% pull(y) %>% na.omit()
wilcox.test(a,b) # p > 0.05
c(mean(a), sd(a), mean(b), sd(b)) %>% round(2)


a <- data_ocean %>% filter(classification == '1') %>% pull(x) %>% na.omit()
b <- data_ocean %>% filter(classification == '3') %>% pull(x) %>% na.omit()
wilcox.test(a,b) # p > 0.05
c(mean(a), sd(a), mean(b), sd(b))


a <- data_ocean %>% filter(classification == '1') %>% pull(length) %>% na.omit() %>% unique()
b <- data_ocean %>% filter(classification == '3') %>% pull(length) %>% na.omit() %>% unique()
t.test(a,b) # p > 0.05
c(mean(a), sd(a), mean(b), sd(b))


data_ocean %>% 
  ggplot(aes(x, y)) + 
  geom_map(map_id = "Australia", map = map, fill = 'grey70') +
  geom_path(data = b2000, colour = 'grey70') +
  geom_point(aes(colour = classification)) + 
  scale_colour_manual(values = c('#ffa600', '#ce5291','#114d78')) +
  labs(x = 'Lon (degrees)', y = 'Lat (degrees)') + 
  lims(x = c(min(data_ocean$x), max(data_ocean$x))) + 
  theme(legend.position = 'right',
                   text = element_text(size = 10),
                   panel.background = element_rect(fill = "white", colour = 'black'),
                   panel.border = element_rect(linetype = "solid", fill = NA))


data_ocean %>% 
  filter(classification != '2') %>% 
  group_by(classification) %>% 
  summarise(unique_seal = length(unique(seal)))




# SIA vs Female SIZE ------------------------------------------------------
size <- read_csv('./data/lnfs_deployment_data.csv')
size <- size %>% dplyr::select(animal_id_1, deploy_mass_kg, deploy_length_cm, deploy_girth_cm) %>% 
  rename(seal = animal_id_1, mass = deploy_mass_kg, length = deploy_length_cm, girth = deploy_girth_cm)

## Fix seal names
i <- which(str_count(size$seal) == 2)
size$seal[i] <- paste0('0', size$seal[i])

## Match size to seal SIA dataframe
h2 <- left_join(h2, size)
h3 <- h2 %>% 
  group_by(seal, year) %>% 
  summarise_at(vars(c13, n15, mass, length, girth), funs(mean(., na.rm  = T)))

h3 %>% 
  ggplot(aes(length, n15)) + 
  geom_point() + 
  geom_smooth() + 
  facet_wrap(~year)

h3 %>% 
  ggplot(aes(length, c13)) + 
  geom_point() + 
  geom_smooth() + 
  facet_wrap(~year)

library(nlme)
library(effects)
library(broom)
lm1 <- gls(n15 ~ length , data = h3)
plot(lm1)
summary(lm1)
par(mfrow=c(2,2))

par(mfrow=c(1,1))

R <- resid(lm1)
F <- fitted(lm1)
lm(R~F) %>% plot()

lm2 <- gls(c13 ~ length  , data = h3, weights = varExp())
summary(lm2)
par(mfrow=c(2,2))
plot(lm2)
par(mfrow=c(1,1))
R <- resid(lm2)
hist(R)
ggqqplot(R)
F <- fitted(lm1)

allEffects(lm1) %>% plot(main = NA, ylab = expression(paste( delta ^{15}, "N ")), xlab = "Length (cm)" )
mtext(expression(paste("(a) ", delta ^{15}, "N")), adj = 0, line = -2, at = c(-45, 0.5))



# Timespent in habitat vs SIZE -------------------------------------------------------
size <- read_csv('./data/lnfs_deployment_data.csv')
size <- size %>% dplyr::select(animal_id_1, deploy_mass_kg, deploy_length_cm, deploy_girth_cm) %>% 
  rename(seal = animal_id_1, mass = deploy_mass_kg, length = deploy_length_cm, girth = deploy_girth_cm)

## Fix seal names
i <- which(str_count(size$seal) == 2)
size$seal[i] <- paste0('0', size$seal[i])

## Match size to seal SIA dataframe
tmp <- left_join(h1, size)

## Add bathy data to get region classification

tmp$z <- raster::extract(bathy_ras, data.frame(x = tmp$x, y = tmp$y))
tmp <- tmp %>% mutate(region = ifelse(z < -2000, 'oceanic', 'shelf'))

tmp <- tmp %>% 
  group_by(seal, region) %>% 
  summarise_at(vars(t.timespent, mass, length, girth), funs(mean(., na.rm  = T)))

lm1 <- gls(t.timespent ~ length *region, data = tmp)
plot(lm1)
summary(lm1)

## Using MCLUST classification instead of bathy
tmp <- data %>% 
  left_join(size) %>% 
  filter(!is.na(t.timespent)) %>% 
  group_by(seal, classification) %>% 
  summarise_at(vars(t.timespent, mass, length, girth), funs(mean(., na.rm  = T)))

lm1 <- gls(log(t.timespent) ~ length *classification, data = tmp)
plot(lm1)
summary(lm1)
plot(allEffects(lm1))

data %>% 
  filter(!is.na(trip)) %>% 
  group_by(seal, trip) %>% 
  summarise(unique_class = length(unique(classification))) %>% 
  pull(unique_class) %>% 
  table # some trips have 2 classification - does that mean seal may have consumed different types of prey in the trip?

# Fit again this time removing trips with 2 classes. 
tmp <- data %>% 
  filter(!is.na(trip)) %>% 
  group_by(seal, trip) %>% 
  mutate(unique_class = length(unique(classification))) %>% 
  filter(unique_class == 1) %>% 
  dplyr::select(seal, trip, classification, year)

tmp <- tmp[!duplicated(tmp),]

tmp <- tmp %>% 
  left_join(h1) %>% 
  left_join(size) %>% 
  group_by(seal, classification) %>% 
  summarise(t.timespent = mean(t.timespent, na.rm =T),
               length = mean(length, na.rm= T), 
               year = first(year))


lm1 <- gls(log(t.timespent) ~ length *classification, data = tmp)
plot(lm1)
summary(lm1)
plot(allEffects(lm1))



# SI between years -------------------------------------------------------

h2 %>% 
  group_by(year) %>% 
  summarise(c13_sd = sd(c13, na.rm = T), n15_sd = sd(n15, na.rm = T),
            c13 = mean(c13, na.rm = T), n15 = mean(n15, na.rm = T) )

h2 %>% 
  group_by(year) %>% 
  filter(root == FALSE) %>% 
  summarise(c13_sd = sd(c13, na.rm = T), n15_sd = sd(n15, na.rm = T),
            c13 = mean(c13, na.rm = T), n15 = mean(n15, na.rm = T) )

## T-test
# Check for normality
shapiro.test(h2$n15)
shapiro.test(h2$c13)

# Non-normal t-test
a <- h2 %>% filter(year == '2016') %>% pull(c13) %>% na.omit()
b <- h2 %>% filter(year == '2017') %>% pull(c13) %>% na.omit()
wilcox.test(a,b) # p <0.05
c(mean(a), sd(a), mean(b), sd(b))


a <- h2 %>% filter(year == '2016') %>% pull(n15) %>% na.omit()
b <- h2 %>% filter(year == '2017') %>% pull(n15) %>% na.omit()
wilcox.test(a,b) # p > 0.05
c(mean(a), sd(a), mean(b), sd(b))




# Non-normal t-test - without root 
a <- h2 %>% filter(root == FALSE) %>% filter(year == '2016') %>% pull(c13) %>% na.omit()
b <- h2 %>% filter(root == FALSE) %>% filter(year == '2017') %>% pull(c13) %>% na.omit()
wilcox.test(a,b) # p <0.05
c(mean(a), sd(a), mean(b), sd(b))

a <- h2 %>% filter(root == FALSE) %>% filter(year == '2016') %>% pull(n15) %>% na.omit()
b <- h2 %>% filter(root == FALSE) %>% filter(year == '2017') %>% pull(n15) %>% na.omit()
wilcox.test(a,b) # p > 0.05
c(mean(a), sd(a), mean(b), sd(b))




# SIA x GLS Temp ----------------------------------------------------------
load("~/Dropbox/phd/LNFS research/LNFS processed data/dive/output/glsTracksTemp.RData")

tr2 <- tr %>% 
  mutate(month = month(Date, abbr = T, label = T)) %>% 
  group_by(id, month) %>% 
  summarise(minSST = mean(MinSST, na.rm = T), minSST_sd = sd(MinSST, na.rm = T))

d3 <- d2 %>% 
  mutate(month = month(date.est, abbr = T, label = T)) %>% 
  group_by(seal, month) %>% 
  summarise(n15 = mean(n15, na.rm = T), c13 = mean(c13, na.rm = T)) %>% 
  rename(id = seal) 

data <- left_join(d3, tr2)
data %>% ggplot(aes(minSST, c13)) +
  geom_point() + 
  geom_smooth() + 
  facet_wrap(~id)

data <- na.omit(data)
hist(data$n15)
hist(data$c13)


library(lme4)
library(nlme)
library(effects)
m1 <- lme(n15 ~ minSST , random = ~ 1 | id, data = data, weights = varExp(form = ~fitted(.)|id), control = lmc)
m2 <- lme(c13 ~ minSST, random = ~ 1 | id, data = data, control = lmc)
m3 <- lme(c13 ~ minSST, random = ~ 1 | id, data = data, control = lmc, weights = varExp())
AIC(m2,m3)
summary(m2)


# SI Plots -------------------------------------------------------------------


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
  mutate(start.date = as.Date(format(start.date, format = '2017-%m-%d')), end.date = as.Date(format(end.date, format = '2017-%m-%d'))) %>% 
  group_by(seal) %>% 
  filter(num.from.base != min(num.from.base))
  
# import gls tracks for comparison
load("./data/glsTracksAll_cleaned.RData")
colnames(tracks) <- tolower(colnames(tracks))
tracks <- tracks %>% mutate(id = as.numeric(as.character(id)))

# plot c13 and n15 time series
# pdf(file = './plots/SI_plots.pdf', width = 10, height = 10)
pal <- wesanderson::wes_palette('BottleRocket2', 3, type = 'continuous' )
t1 <- min(d3$start.date) 
t2 <- max(d3$end.date)

mytheme <- theme(legend.position = 'bottom',
                 text = element_text(size = 10),
                 panel.background = element_rect(fill = "white", colour = 'black'),
                 panel.border = element_rect(linetype = "solid", fill = NA))
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
  geom_line(size = 0.5) +
  geom_point(size = 1, aes(colour = has.glycine)) + 
  geom_vline(aes(xintercept = start.date), linetype = 'dashed') +
  geom_vline(aes(xintercept = end.date), linetype = 'dashed') +
  labs(x = 'Date', y = expression(paste( delta ^"13", "C"))) +
  scale_x_date(date_breaks = '1 month', date_labels = '%b', limits = c(t1, t2)) +
  mytheme + 
  facet_wrap(~seal) + 
  scale_colour_manual(values = wesanderson::wes_palette('Rushmore1')[c(3,5)], name = 'Glycine?')

# date vs d15n (without glycine segments)
d3 %>% 
  mutate(n15_2 = ifelse(has.glycine == TRUE, NA, n15 )) %>% 
  ggplot(aes(x = date, y = n15_2)) + 
  geom_point(size = 1) + 
  geom_line(size = 0.5) +
  geom_vline(aes(xintercept = start.date), linetype = 'dashed') +
  geom_vline(aes(xintercept = end.date), linetype = 'dashed') +
  labs(x = 'Date', y = expression(paste( delta ^"15", "N"))) +
  scale_x_date(date_breaks = '1 month', date_labels = '%b', limits = c(t1,t2)) +
  mytheme + 
  facet_wrap(~seal)

# date vs d15n (only glycine seals)
d3 %>% 
  filter(glycine == TRUE) %>% 
  ggplot(aes(x = len.from.base, y = n15)) +
  geom_line(size = 0.5) +
  geom_point(size = 1, aes(colour = has.glycine)) + 
  # geom_vline(aes(xintercept = start.date)) +
  # geom_vline(aes(xintercept = end.date)) +
    geom_hline(aes(yintercept = 17.6), linetype = 'dashed') +
  labs(x = 'Length from base', y = expression(paste( delta ^"15", "N")) ) +
  # scale_x_date(date_breaks = '1 month', date_labels = '%b') +
  mytheme +  
  facet_wrap(~seal) + 
  scale_colour_manual(values = wesanderson::wes_palette('Rushmore1')[c(3,5)], name = 'Glycine?')
# dev.off()



pdf(file = './plots/individual_SI_plots_autoMethod.pdf', width = 10, height = 7)
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
  geom_boxplot() 
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
