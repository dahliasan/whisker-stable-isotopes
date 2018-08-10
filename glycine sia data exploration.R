# explore stable isotope results

library(tidyverse)
library(zoo)
library(wesanderson)


# read in data
d <- read_csv("./glycine animals/lnfs_glycine_whisker_bulk_sia_correctedsequence.csv") %>% 
  select(seal, glycine.start.date, date, num, root, n15, c13, sample.len)

d1 <- d %>% mutate(missing = ifelse(is.na(n15), TRUE, FALSE))


d1 %>% ggplot(aes(x = num, group = missing)) +
  geom_line(aes(y = n15, colour = root)) + 
  geom_point(aes(y = n15), size = 0.5) +
  facet_wrap(~seal)

# natural range of n15 for lactating lnfs at cape gantheaume is max 16. therefore glycine spikes are values > 16
# assign glycine injection dates to the first sample whose n15 > 17
d1 <- d1 %>% mutate(glycine = ifelse(n15 > 17.5, TRUE, FALSE)) 

# interpolate dates for samples in between glycine injection and whisker collection dates

# d1 <- d1 %>% mutate(date = as.Date(ifelse(glycine == TRUE, glycine.start.date, date))) %>% 
#   arrange(seal, num)

g <- d1 %>% 
  split(.$seal) %>% 
  map(~.$num[max(which(.$n15 > 17.5))])

map2(split(d1, d1$seal), g, function(x, y){
  x$date[y] <- x$glycine.start.date[1]
  return(x)
}) %>% reduce(bind_rows) -> d1

d1 %>% arrange(seal, num) %>% 
  split(.$seal) %>% 
  map(function(x){
    t <- x$date
    tt <- as.Date(na.approx(t))
    x$date[1:length(tt)] <- tt
    return(x)
  }) %>% reduce(bind_rows) -> d1

d1$seal <- as.factor(d1$seal)

# save plots 
pdf(file = 'lnfs_glycine_SIA.pdf')
d1 %>% mutate(num = num * -1) %>% 
  ggplot(aes(x = num, group = missing)) +
  geom_line(aes(y = n15)) + 
  geom_point(aes(y = n15, colour = glycine), size = 1) +
  facet_wrap(~seal) +
  scale_colour_manual(values = wes_palette('Darjeeling1', 2))

d1 %>% mutate(num = num * -1) %>% 
  ggplot(aes(x = num, group = missing)) +
  geom_line(aes(y = c13)) + 
  geom_point(aes(y = c13, colour = glycine), size = 1) +
  facet_wrap(~seal) +
  scale_colour_manual(values = wes_palette('Darjeeling1', 2))

d1 %>% filter(n15 <= 17.5, missing == FALSE) %>% 
  ggplot(aes(date, n15, group = missing)) + 
  geom_line() +
  geom_point( size = 1) + 
  facet_wrap(~seal) + 
  scale_x_date(date_breaks = '1 month', date_labels = '%b')

d1 %>% 
  ggplot(aes(date, c13)) + 
  geom_line() +
  geom_point( size = 1) + 
  facet_wrap(~seal) + 
  scale_x_date(date_breaks = '1 month', date_labels = '%b')

dev.off()

# calculate growth rate
d1 %>% group_by(seal) %>% 
  mutate(difftime =  lag(date) - lead(date)) %>% 
  summarise(growth.rate = round(mean(difftime, na.rm = T)/2))

