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
# library(ggthemr)
library(car)
library(MASS) 
library(sjPlot)
library(lme4)
library(nlme)
library(MuMIn)
library(DHARMa)
library(effects)
library(broom)
library(mclust)
library(ellipse)
options(dplyr.width = Inf) #enables head() to display all coloums

## Load data
# adults analysed were regrowths
# pup whiskers were plucked at the end of fieldwork
d0 <- read_csv("./data/lnfs_all_whisker_bulk_sia_correctedsequence.csv")
d1 <- d0 %>% 
  filter(adult == TRUE) %>% 
  dplyr::select(seal, num.from.base, c13, n15, sample.len, root, root.length, original.length, plucked, start.date, end.date, glycine) %>% 
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
  dplyr::select(-norootlen)

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
  dplyr::select(seal, num.from.base, c13, n15, sample.len, root, root.length, original.length, plucked, start.date, end.date, glycine, growth.rate, res) %>% 
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
  dplyr::select(seal, gly.start)
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
load("C:/Users/footoddy/Documents/phd/LNFS research/LNFS processed data/animal home range/time spent/timespent_and_percentage_by_seal_season_trip.RData")

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

# Extract TOPO (region)  ---------------------------------------------
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

## Assign foraging state
h2$state <- 'provisioning'
h2$state[h2$seal %in% c('307', '317', '324', '340')] <- 'unconstrained'
h2$state <- as.factor(h2$state)

## Check if region assigned correctly: 
## Should mostly have 069 and 071 seals without region because they don't have tracks. 
h2 %>% filter(is.na(region)) %>% dplyr::select(seal, region) %>% View()

h2 %>% filter(!is.na(region)) %>% 
  ggplot(aes(x = region, y = n15,  group = region)) + geom_boxplot()
h2 %>% filter(!is.na(region)) %>% 
  ggplot(aes(x = region, y = c13,  group = region)) + geom_boxplot()

# Load female SIZE data ----RUN CODE UP TO HERE.-----------------------------------------------
size <- read_csv('./data/lnfs_deployment_data.csv')
size <- size %>% dplyr::select(animal_id_1, deploy_mass_kg, deploy_length_cm, deploy_girth_cm) %>% 
  rename(seal = animal_id_1, mass = deploy_mass_kg, length = deploy_length_cm, girth = deploy_girth_cm)

## Fix seal names
i <- which(str_count(size$seal) == 2)
size$seal[i] <- paste0('0', size$seal[i])

# Latitudal pattern -------------------------------------------------------

mytheme <- theme(legend.position = 'bottom',
                 text = element_text(size = 8.5),
                 panel.background = element_rect(fill = "white", colour = 'black'),
                 panel.border = element_rect(linetype = "solid", fill = NA))

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
plot(m2$gam, select = 1, xlab = NA, ylab = 's(lat,1):year2016', scheme = 1 )
mtext(expression(paste("(a) ", delta ^"13", "C")), adj = 0, line = -2, at = c(-45, 0.5))
plot(m2$gam, select = 2, xlab = NA, ylab = 's(lat,1):year2017', scheme = 1 )
mtext(expression(paste("(b) ", delta ^"13", "C")), adj = 0, line = -2, at = c(-45, 0.5))
plot(m1$gam, select = 2, xlab = 'Latitude (degrees)', ylab = 's(lat,3.58):year2017', scheme = 1)
mtext(expression(paste("(c) ", delta ^"15", "N")), adj = 0, line = -2, at = c(-45, 0.5))

# SI vs Region + Year ----------------------------------------------------
# Fit GLM
library(lme4)
library(broom)
library(effects)
library(bbmle)
library(finalfit)

## -------- Including root sections in 2017 
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
AICc(m0,m1,m2) 
summary(m2)
plot(m2)
hist(residuals(m2))

ICtab(m0,m1,m2, base = T, logLik = T, type = 'AIC') 
m2 %>% fit2df()

## 15N

mm0 <- glm(n15 ~ region*year, data = glm_data)
mm1<- lmer(n15 ~ region*year + (1|seal), data = glm_data)
anova(mm1, mm0) # yes keep seal random effect
drop1(mm1, test="Chisq")
mm2 <- update(mm1, .~. - region:year)
drop1(mm2, test = 'Chisq')
ICtab(mm0,mm1,mm2, base = T, logLik = T, type = 'AIC') 
summary(mm1)
plot(allEffects(mm1))

# write.csv(rbind(tidy(m2), tidy(mm1)), './output/GLM_final_fit.csv')

## -------- Excluding root sections in 2017
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

# write.csv(rbind(tidy(m2), tidy(mm1)), './output/GLM_final_fit_exclude_root.csv')

# Effect plots
library(sjPlot)
library(ggpubr)
mytheme <- theme(
  text = element_text(size = 11),
  panel.background = element_rect(fill = "white", colour = 'black'),
  panel.border = element_rect(linetype = "solid", fill = NA)
)
p1 <- plot_model(m2, type = 'eff', terms = 'region', show.ci = TRUE) + 
  lims(y = c(-16.1 , -15.65)) + 
  mytheme + 
  labs(x = 'Habitat', y = expression(paste( delta ^"13", "C")), title = NULL)
p2 <- plot_model(m2, type = 'eff', terms = 'year', show.ci = TRUE) + 
  lims(y = c(-16.1 , -15.65)) + 
  mytheme + 
  labs(x = 'Year', y = expression(paste( delta ^"13", "C")), title = NULL)
pp1 <- ggarrange(p1, p2, ncol = 2, labels = c('a', 'b'))

p3 <- plot_model(mm1, type = 'int', grid = T, colors = c('black', 'black')) + 
  mytheme + 
  labs(x = 'Habitat', y = expression(paste( delta ^"15", "N")), title = NULL)

# pdf(file = './plots/SI_vs_region_year.pdf', height = 6, width = 6)
ggarrange(pp1, p3, nrow = 2, labels = c('', 'c'))

# dev.off()

## Only provisioning females displaying shelf + oceanic foraging + exclude root
glm_data <- h2 %>% 
  filter(!is.na(region), root == FALSE) %>% 
  filter(!seal %in% c('329', '353', '315', '318') ) %>% 
  filter(state == 'provisioning') %>% 
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
AICc(mm0,mm1,mm2)
summary(mm1)
plot(allEffects(mm1))


plot_model(mm1, type = 'int', grid = T, colors = c('black', 'black')) + 
  mytheme + 
  labs(x = 'Habitat', y = expression(paste( delta ^"15", "N")), title = NULL)


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

# Assign foraging state
data <- left_join(data, h2) 
data$fs <- 'shelf/oceanic'
data$fs[data$seal %in% c('329', '353')] <- 'oceanic'
data$fs[data$seal %in% c('318', '315')] <- 'shelf'
data$fs[data$seal %in% c('307', '324', '317', '340')] <- 'unconstrained'
data <- data %>% mutate(year = as.factor(year(start.date)))
data$uncertainty <- m1$uncertainty

# ------------ Biplot
data$region[is.na(data$region)] <- "NA"
p1 <- data %>%
  ggplot(aes(c13, n15)) +
  # geom_point(aes(shape = classification, colour = region, size = uncertainty))+
  geom_point(aes(shape = region, colour = classification, size = uncertainty))+
  geom_path(data = mclust.el, aes(group = classification, linetype = classification)) + 
  theme_pubr(base_size = 11, legend = "right") +  
  # scale_color_brewer(type = 'qual', palette = 7, direction = -1) + 
  scale_colour_manual(values = c('#ffa600', '#ce5291', '#114d78', 'red', 'blue', 'green')) + 
  scale_shape_manual(values = c(1, 16, 3)) + 
  labs(x = expression(paste(delta ^"13", "C")), y = expression(paste(delta ^"15", "N")))

# Prey biplot --------------------------------------------------------------  
# load nearby prey/organism SI
prey <- read_csv('./data/fish_si.csv') %>% 
  filter(lnfs.related == TRUE) %>% 
  mutate(n15_se = n15_sd/sqrt(n), c13_se = c13_sd/sqrt(n))

# pdf('./plots/biplot_mclust_prey.pdf', width = 7.5, height = 7)
p1 + 
  geom_point(data = prey, aes(c13, n15, colour = habitat) ) + 
  geom_text(data = prey, aes(label = shortname), colour = 'black', size = 3, nudge_y = -0.2) +
  geom_errorbar(data = prey, aes(x = c13, xmin = c13-c13_se, xmax = c13+c13_se, 
                                 y = n15, ymin = n15-n15_se, ymax = n15+n15_se), width = 0.05)
# dev.off()

# ------- Map of classication
library(marmap)
map <- map_data("worldHires", "Australia")
bathy <- marmap::getNOAA.bathy(129,151,-46, -35,res=1, keep=TRUE)
bathy_map <- tbl_df(fortify(bathy))
b2000 <- bathy_map %>% filter(z >= -2000 & z <= -1000) %>% group_by(x) %>% summarise(y = mean(y))

mytheme <- theme(legend.position = 'right',
                 text = element_text(size = 10),
                 panel.background = element_rect(fill = "white", colour = 'black'),
                 panel.border = element_rect(linetype = "solid", fill = NA))
data_map <- data %>% 
  mutate(yday = yday(date.est)) %>% 
  filter(!is.na(x)) %>%
  left_join(size) %>% 
  dplyr::select(seal, yday, classification, year, length, x, y)


p2 <- data_map %>% 
  ggplot(aes(x, y)) + 
  geom_map(map_id = "Australia", map = map, fill = 'grey70') +
  geom_path(data = b2000, colour = 'grey70') +
  geom_point(aes(colour = classification), size = 2) +
  # scale_colour_viridis_c() + 
  scale_colour_manual(values = c('#ffa600', '#ce5291','#114d78')) +
  labs(x = 'Lon (degrees)', y = 'Lat (degrees)') + 
  lims(x = c(min(data_map$x), max(data_map$x))) + 
  facet_wrap(~year) + 
  theme(legend.position =  c(0.12, 0.2),
        text = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = 'black'),
        panel.border = element_rect(linetype = "solid", fill = NA))

# pdf('./plots/mclust_biplot_spatialmap.pdf', width = 11, height = 4.5)
ggarrange(p1, p2, widths = c(1, 1.1), labels = c('a', 'b'))
# dev.off()

## ---------- Uncertainty
q <- quantile(m1$uncertainty)
uncertainty <- tbl_df(data.frame(uncertainty = m1$uncertainty, quantile = NA))
uncertainty$quantile <- cut(uncertainty$uncertainty, q, labels = c('q25', 'q50', 'q75', 'q100'))
uncertainty$quantile[which(is.na(uncertainty$quantile))] <- 'q25'
m1$uncertainty %>% hist()
plot(m1, what = 'uncertainty')

## --------------- Classification by year
data <- data %>%
  mutate(state = ifelse(fs == 'unconstrained', 'unconstrained', 'provisioning'))
table(data$classification, data$year, data$state)

data %>%
  ggplot(aes(yday, classification)) + 
  geom_point() + 
  facet_wrap(~year)


# Individual fidelity -----------------------------------------------------

# pdf('./plots/individual_isotope_signatures_clusters.pdf', width = 7, height = 7)
data %>% 
  ggplot(aes(c13, n15)) + 
  geom_point(aes(colour = classification)) +
  scale_colour_manual(values = c('#ffa600', '#ce5291','#114d78')) +
  facet_wrap(year~state+seal) + 
  labs(x = expression(paste(delta ^"13", "C")), y = expression(paste(delta ^"15", "N"))) +
  mytheme

# dev.off()
# Cluster vs TIme ---------------------------------------------------------

## Classification vs month (aka season)
data$month <- month(data$date.est)
data$yday <- yday(data$date.est)
data_ocean <- data %>% 
  filter(classification != '2') %>% 
  mutate(cluster = ifelse(classification == '3', '1', '0'), cluster = as.factor(cluster))

data_ocean$classification <- droplevels(data_ocean$classification)

data_ocean %>% 
  group_by(classification) %>% 
  summarise(mean(yday))

library(DHARMa)
m0 <- glm(cluster~ log(yday)*year, family = binomial, data = data_ocean)
m1 <- glmer(cluster ~ log(yday)*year + (1|seal) , family = binomial, data = data_ocean)
anova(m1, m0)
drop1(m1, test = 'Chisq' )
simulationOutput <- simulateResiduals(fittedModel = m1)
plotSimulatedResiduals(simulationOutput = simulationOutput) # you want p value > 0.05
# test uniform residuals you want p value > 0.05, it also gives some indication of overdispersion (if p < 0.05)
testUniformity(simulationOutput = simulationOutput) 
plot(m1)
summary(m1)

# pdf(file = './plots/ocean_classification_vs_yday_year.pdf', width = 5, height = 4)
plot(allEffects(m1), main = '', ylab = 'Probability of cluster 3', xlab = 'Day of the year')
# dev.off()

## only provisioning females
data_ocean_subset <- data_ocean %>% 
  filter(fs !='unconstrained')
m1 <- glmer(classification ~ log(yday) + (1|seal) , family = binomial, data = data_ocean_subset)
drop1(m1, test = 'Chisq' )
simulationOutput <- simulateResiduals(fittedModel = m1)
plotSimulatedResiduals(simulationOutput = simulationOutput) # you want p value > 0.05
# test uniform residuals you want p value > 0.05, it also gives some indication of overdispersion (if p < 0.05)
testUniformity(simulationOutput = simulationOutput) 
plot(m1)
summary(m1)


# n clusters per female ---------------------------------------------------
## How many different clusters for each female?
data %>% 
  group_by(seal) %>% 
  summarise(length(unique(classification))) 

## Plot timeline of cluster classification for each individual

# This stops the line from joining completely for 351 cause there are missing datapoints.
u <- data %>% filter(seal == '351') %>% slice(1)
u <- u %>% mutate(trip = '2', n15 = NA, c13 = NA)
tmp <- bind_rows(data, u) %>% 
  arrange(seal, yday)
  
# plot
# pdf('./plots/Cluster_classification_timeline_individual.pdf', width = 10, height = 11)
tmp %>% 
  mutate(state = ifelse(seal %in% c('307', '317', '324', '340'), 'unconstrained', 'provisioning')) %>% 
  ggplot(aes(yday)) + 
  geom_line(aes(y = n15)) + 
  geom_point(aes(y = n15, colour = classification, size = uncertainty)) +
  geom_line(aes(y = c13 * -1), linetype = 'dashed') + 
  geom_point(aes(y = c13 * -1, colour = classification, size = uncertainty), shape = 18) +
  scale_colour_manual(values = c('#ffa600', '#ce5291','#114d78'), name = 'cluster') +
  scale_y_continuous(sec.axis = sec_axis(~.*-1, name = expression(paste( delta ^{13}, "C ")))) + 
  facet_wrap(state~fs+seal ) + 
  theme_bw() + 
  labs(x = 'Day of year', y = expression(paste( delta ^{15}, "N ")))
# dev.off()

## Uncertainty of classification against time?

data %>% 
  ggplot(aes(yday, uncertainty)) + 
  geom_smooth()

m1 <- lmer(uncertainty ~ yday + (1|seal), data = data)
drop1(m1, test = 'Chisq') # No relationship 



# SIA vs Female SIZE ------------------------------------------------------
## Match size to seal SIA dataframe
h3 <- left_join(h2, size)
h3 <- h3 %>% 
  group_by(seal, year, state) %>% 
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
lm1 <- gls(n15 ~ length*state , data = h3)

plot(lm1)
summary(lm1)

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

## Only provisioning females
h3.subset <- h3 %>% 
  filter(state == 'provisioning')
lm1 <- gls(n15 ~ length , data = h3.subset)
plot(lm1)
summary(lm1)

lm1 <- gls(c13 ~ length , data = h3.subset)
plot(lm1)
summary(lm1)



# Timespent in habitat vs size / year -------------------------------------------------------
library(MASS)
library(MuMIn)

## Match size to seal SIA dataframe
tmp <- left_join(h1, size)

## Add bathy data to get region classification

tmp$z <- raster::extract(bathy_ras, data.frame(x = tmp$x, y = tmp$y))
tmp <- tmp %>% mutate(region = ifelse(z < -2000, 'oceanic', 'shelf'))

tmp <- tmp %>%
  # filter(!seal %in% c('307', '317', '324', '340')) %>% 
  group_by(seal, region) %>%
  mutate(year = year(median.date)) %>%
  summarise_at(vars(t.timespent, mass, length, girth, year), funs(mean(., na.rm  = T)))

tmp$state <- 'provisioning'
tmp$state[tmp$seal %in% c('307', '317', '324', '340')] <- 'unconstrained'
tmp$state <- as.factor(tmp$state)

tmp$year <- as.factor(tmp$year)
tmp$region <- as.factor(tmp$region)

tmp <- tmp %>% filter(state == 'provisioning', !seal %in% c('353','329','
                                                            315','318'))

lm1 <- lm(sqrt(t.timespent) ~ length * region + region*year, data = tmp)
lm2 <- lm(sqrt(t.timespent) ~ length + region + region*year , data = tmp)
lm3 <- lm(sqrt(t.timespent) ~ length + region + year , data = tmp)
lm4 <- lm(sqrt(t.timespent) ~ region + length , data = tmp)
lm5 <- lm(sqrt(t.timespent) ~ region  , data = tmp)
AICc(lm1, lm2, lm3, lm4, lm5) 

# Diagnostics
library(car)
outlierTest(lm4)
studres(lm4) %>% hist()
qqPlot(lm4)
ncvTest(lm4)
spreadLevelPlot(lm4)
plot(lm4)
summary(lm4)

plot(allEffects(lm4))

p1 <- plot_model(lm2, type = 'eff', colors = 'bw', show.loess = F, term = 'region') +
  mytheme + 
  labs(y = 'Log(cell residence time) (h)', x = 'Habitat', title = NULL) 
p2 <- plot_model(lm2, type = 'eff', colors = 'bw', show.loess = F, term = 'length') +
  mytheme + 
  labs(y = NULL, x = 'Length (cm)', title = NULL) 

# pdf(file = './plots/timespent_habitat_length.pdf', width = 6.5, height = 3)
p3 <- ggarrange(p1, p2, labels = c('b', 'c'))
# dev.off()

## ------- only shelf + oceanic provisioning females 
tmp <- tmp %>% filter(state == 'provisioning', !seal %in% c('353','329','
                                                            315','318'))


lm1 <- lm(sqrt(t.timespent) ~ length * region + region*year, data = tmp)
summary(lm1)
lm2 <- lm(sqrt(t.timespent) ~ length + region + region*year , data = tmp)
summary(lm2)
lm3 <- lm(sqrt(t.timespent) ~ region * year , data = tmp)
summary(lm3)
lm4 <- lm(sqrt(t.timespent) ~ region + year , data = tmp)
summary(lm4)
lm5 <- lm(sqrt(t.timespent) ~ region  , data = tmp)
AICc(lm1, lm2, lm3, lm4, lm5) 



## ------- Using MCLUST classification instead of bathy ##
# tmp <- data %>% 
#   left_join(size) %>% 
#   filter(!is.na(t.timespent)) %>% 
#   group_by(seal, classification) %>% 
#   mutate(year = year(median.date)) %>% 
#   summarise_at(vars(t.timespent, mass, length, girth, year), funs(mean(., na.rm  = T)))
# 
# lm1 <- lm(t.timespent ~ length *classification, data = tmp)
# lm2 <- stepAIC(lm1, direction = 'backward')
# lm2$anova
# plot(lm2)
# summary(lm2)
# View(tidy(lm2))
# # finalfit_list[[2]] <- tidy(lm2)
# plot(allEffects(lm2))

data %>% 
  filter(!is.na(trip)) %>% 
  group_by(seal, trip) %>% 
  summarise(unique_class = length(unique(classification))) %>% 
  pull(unique_class) %>% 
  table # some trips have 2 classification - does that mean seal may have consumed different types of prey in the trip?


data %>% 
  filter(!is.na(trip)) %>% 
  group_by(seal, trip) %>% 
  summarise(unique_class = length(unique(classification))) %>% 
  filter(unique_class >1)

## ------- Cluster groups: only trips with 1 class ##
tmp <- data %>% 
  filter(!is.na(trip)) %>% 
  group_by(seal, trip) %>% 
  mutate(unique_class = length(unique(classification))) %>% 
  filter(unique_class == 1) %>% 
  dplyr::select(seal, trip, classification, year, fs) %>% 
  mutate(state = ifelse(fs == 'unconstrained', 'unconstrained', 'provisioning')) %>% 
  ungroup() %>% 
  mutate(state = as.factor(state)) %>% 
  filter(classification != '2') # remove cluster 2 since it has low number of samples.

tmp <- tmp[!duplicated(tmp),]

tmp <- tmp %>% 
  left_join(h1) %>% 
  left_join(size) %>% 
  group_by(seal, classification, state) %>% 
  summarise(t.timespent = mean(t.timespent, na.rm =T),
            length = mean(length, na.rm= T), 
            year = first(year))

lm1 <- lm(sqrt(t.timespent) ~ length * classification  ,data = tmp)
lm2 <- lm(sqrt(t.timespent) ~ length + classification ,data = tmp)
AICc(lm1, lm2)

summary(lm1)
# View(tidy(lm1))

# Diagnostics
library(car)
outlierTest(lm1)
studres(lm1) %>% hist()
qqPlot(lm1)
ncvTest(lm1)
spreadLevelPlot(lm1)
# plot(lm2)

# pdf(file = './plots/timespent_oceanClassification_length_1classTripsOnly.pdf', width = 4, height = 3)
plot_model(lm1, type = 'int', grid = T, color = 'bw') + 
  mytheme + 
  labs(y = 'Sqrt(cell residence time) (h)', x = 'Length (cm)', title = NULL) 
# dev.off()


## ---------  only provisioning females ##
tmp$state <- 'provisioning'
tmp$state[tmp$seal %in% c('307', '317', '324', '340')] <- 'unconstrained'
tmp$state <- as.factor(tmp$state)

tmp <- tmp %>% 
  filter(state == 'provisioning')
tmp_outrm <- tmp[-10, ]

lm1 <- lm(sqrt(t.timespent) ~ length *classification,data = tmp_outrm)
lm2 <- lm(sqrt(t.timespent) ~ length + classification,data = tmp_outrm)
lm3 <- lm(sqrt(t.timespent) ~ classification, data = tmp_outrm)
AICc(lm1, lm2,lm3)
summary(lm3)
View(tidy(lm3))

# Diagnostics
library(car)
outlierTest(lm3)
studres(lm3) %>% hist()
qqPlot(lm3)
ncvTest(lm3)
spreadLevelPlot(lm3)
plot(lm3)


finalfit_list[[3]] <- tidy(lm2)
View(tidy(lm2))

# pdf(file = './plots/timespent_classification_length_1classTripsOnly_provisioningOnly.pdf', width = 3, height = 3)
plot_model(lm3, type = 'eff', terms = 'classification', color = 'bw') + 
  mytheme + 
  labs(y = 'Cell residence time (h)', x = 'Cluster', title = NULL) 
# dev.off()



# Summarise for manuscript ------------------------------------------------
summ <- d2 %>% group_by(seal) %>% 
  summarise(startdate = first(start.date),
            enddate = first(end.date),
            glycine = first(glycine), 
            growthrate = first(growth.rate),
            tlength = first(original.length),
            rlength = first(root.length),
            nsections = n())


tmp <- h %>% 
  group_by(seal, trip) %>%
  summarise(tripdur = round(difftime(first(enddate), first(startdate), units = 'days'), 2)) 

left_join(h2, tmp) %>% 
  filter(!is.na(tripdur), region =='oceanic') %>% 
  group_by(seal) %>% 
  summarise(tripdur = mean(tripdur)) %>% 
  left_join(summ) -> tmp
  
lm(tlength~tripdur, data = tmp) %>% summary()


h1 %>% 
  mutate(state = ifelse(seal %in% c('307', '317', '324', '340'), 'unconstrained', 'provisioning'),
         state = as.factor(state)) %>% 
  group_by(seal, state) %>%
  summarise(t.timespent = mean(t.timespent, na.rm = T)) %>% 
  left_join(summ) %>% 
  # filter(state == 'provisioning') %>% 
  ggplot(aes(tlength, t.timespent)) + 
  # geom_boxplot() + 
  geom_smooth() + 
  geom_point(aes(colour =state ))


## Test between glycine vs non-glycine growth rates
a <- summ %>% pull(growthrate) # non-glycine method
b <- c(0.186, 0.271, 0.136, 0.246, 0.186, 0.151, 0.155, 0.17)

wilcox.test(a,b) # p  > 0.05
round(c(mean(a), sd(a), mean(b), sd(b)),2)


## isotope value summary for each individual

tmp <- d2 %>% 
  group_by(seal) %>% 
  filter(has.glycine == FALSE) %>% 
  summarise(mean(n15, na.rm = T), sd(n15, na.rm = T), min(n15, na.rm = T), max(n15, na.rm = T),
            mean(c13, na.rm = T), sd(c13, na.rm = T), min(c13, na.rm = T), max(c13, na.rm = T))

tmp <- d2 %>% 
  filter(has.glycine == FALSE) %>% 
  summarise(mean(n15, na.rm = T), sd(n15, na.rm = T), min(n15, na.rm = T), max(n15, na.rm = T),
            mean(c13, na.rm = T), sd(c13, na.rm = T), min(c13, na.rm = T), max(c13, na.rm = T)) %>% 
  bind_rows(tmp)

tmp <- data %>% 
  group_by(classification) %>% 
  summarise(mean(n15, na.rm = T), sd(n15, na.rm = T), min(n15, na.rm = T), max(n15, na.rm = T),
            mean(c13, na.rm = T), sd(c13, na.rm = T), min(c13, na.rm = T), max(c13, na.rm = T)) %>% 
  bind_rows(tmp)
# write.csv(tmp, './output/individual_isotope_summary.csv')

## Test SI of provsioning vs unconstrained females
# c13
a <- h2 %>% filter(!is.na(c13), state == 'provisioning') %>% pull(c13)
b <- h2 %>% filter(!is.na(c13), state == 'unconstrained') %>% pull(c13)
wilcox.test(a,b) # p  > 0.05
round(c(mean(a), sd(a), mean(b), sd(b)),2)

# n15
a <- h2 %>% filter(!is.na(n15), state == 'provisioning') %>% pull(n15)
b <- h2 %>% filter(!is.na(n15), state == 'unconstrained') %>% pull(n15)
shapiro.test(a)
shapiro.test(b)

t.test(a,b) # p  > 0.05
round(c(mean(a), sd(a), mean(b), sd(b)),2)


## Size of provisioning vs unconstrained females
s <- left_join(h2, size) %>% 
  group_by(seal) %>% 
  summarise(length = first(length)) %>% 
  mutate(state = ifelse(seal %in% c('307', '317', '324', '340'), 'unconstrained', 'provisioning'),
         state = as.factor(state)) 

a <- s %>% filter(state == 'provisioning') %>% pull(length)
b <- s %>% filter(state == 'unconstrained') %>% pull(length)
t.test(a, b, 'two.sided') # t = 2.5494, df = 5.6859, p-value = 0.04565
round(c(mean(a), sd(a), mean(b), sd(b)),2)

## d13C on the shelf only
s <- h2 %>% filter(!is.na(region), root == FALSE) %>% 
  mutate(region = as.factor(region))

a <- s %>% filter(region == 'shelf', year == '2016') %>% pull(c13) %>% na.omit
b <- s %>% filter(region == 'shelf', year == '2017') %>% pull(c13) %>% na.omit

shapiro.test(a); shapiro.test(b)
t.test(a, b, 'two.sided') # t = 2.7632, df = 22.801, p-value = 0.01112
round(c(mean(a,), sd(a), mean(b), sd(b)),2)

# D15N
a <- s %>% filter(region == 'oceanic', year == '2016') %>% pull(c13) %>% na.omit
b <- s %>% filter(region == 'oceanic', year == '2017') %>% pull(c13) %>% na.omit
shapiro.test(a); shapiro.test((b))
wilcox.test(a,b) # p  > 0.05
round(c(mean(a), sd(a), mean(b), sd(b)),2)



# Plot study region map ---------------------------------------------------
## Plot tracks showing core foraging areas
library(marmap)
map <- map_data("worldHires", "Australia")
bathy <- marmap::getNOAA.bathy(129,151,-46, -35,res=1, keep=TRUE)
bathy_map <- tbl_df(fortify(bathy))
b2000 <- bathy_map %>% filter(z >= -2000 & z <= -1000) %>% group_by(x) %>% summarise(y = mean(y))

mytheme <- theme(legend.position = 'right',
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

# Prepare data 
smap.dat <- h2 %>% 
  mutate(month = month(median.date, abbr = T, label = T)) %>% 
  group_by(x, y, year, state) %>% 
  summarise(c13_sd = sd(c13, na.rm = T), n15_sd = sd(n15, na.rm = T),
            c13 = mean(c13, na.rm = T), n15 = mean(n15, na.rm = T))

# Overall mean spatial map for each year

library(jcolors)
p1 <- smap.dat %>%
  ungroup() %>% 
  mutate(x = round(x, 10), y = round(y,10)) %>% 
  filter(!is.na(n15)) %>% 
  ggplot(aes(x, y))  +
  geom_point(aes(alpha = n15, colour = n15), size = 1) + 
  geom_map(map_id = "Australia", map = map, fill = 'grey70') +
  geom_path(data = b2000, colour = 'grey70', size = 0.3) +
  facet_wrap(year~state) + 
  scale_colour_viridis(option = 'A', direction = -1,  name = expression(paste( delta ^"15", "N", ' (‰)'))) +
  scale_alpha_continuous(guide = 'none', name = expression(paste( delta ^"15", "N"))) +
  # scale_colour_jcolors_contin('pal12') +  
  mytheme + 
  lims(x = c(min(smap.dat$x, na.rm = T), max(smap.dat$x, na.rm = T))) + 
  labs(x = NULL, y = 'Lat (degrees)')

p2 <- smap.dat %>%
  ungroup() %>% 
  mutate(x = round(x, 10), y = round(y,10)) %>% 
  filter(!is.na(c13)) %>% 
  ggplot(aes(x, y))  +
  geom_point(aes(alpha = c13, colour = c13), size = 1, alpha = 0.8) + 
  geom_map(map_id = "Australia", map = map, fill = 'grey70') +
  geom_path(data = b2000, colour = 'grey70', size = .3) +
  facet_wrap(year~state) + 
  scale_colour_viridis(option = 'A', direction = -1,  name = expression(paste( delta ^"13", "C", " (‰)"))) +
  scale_alpha_continuous(guide = 'none', name = expression(paste( delta ^"13", "C"))) +
  mytheme +
  lims(x = c(min(smap.dat$x, na.rm = T), max(smap.dat$x, na.rm = T))) + 
  labs(x = 'Lon (degrees)', y = 'Lat (degrees)')

# pdf(file = './plots/SI_spatialmap_year_foragingstate.pdf', width = 6, height = 5)
ggarrange(p1, p2, nrow = 2)
# dev.off()

## -------- Show all core foraging cells instead of just the median location
a <-  h %>% 
  group_by(seal, trip, median.date) %>%
  filter(timespent >= q) %>% 
  filter(y <= min(y) + 2)

b <- h2  %>% 
  dplyr::select(seal, median.date, n15, c13, year, state)

smap.dat <- left_join(a, b) %>% 
  filter(!is.na(state)) %>% 
  group_by(x, y, year, state) %>% 
  summarise(n15 = mean(n15, na.rm = T), c13 = mean(c13, na.rm = T))



smap.dat %>% ungroup() %>% 
  filter(!is.na(n15)) %>% 
  ggplot(aes(x, y))  +
  geom_raster(aes(fill = n15)) + 
  geom_map(map_id = "Australia", map = map, fill = 'grey70') +
  geom_path(data = b2000, colour = 'grey70') +
  facet_wrap(year~state) + 
  scale_fill_viridis(option = 'D') +
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

# dev.off()

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



# pdf(file = './plots/individual_SI_plots_autoMethod.pdf', width = 10, height = 7)
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
