## Estimate dates along regrowth by 1) identifying actual regrowth section of whisker (i.e. ignore under skin part of whisker which will be estimated)
# 2) estimate dates by first identifying the start date and forward estimate (deployment period) and backward (predeploy) estimate

estRegrowthDates <- function(x, ratio = 1, method = "auto", estimateBackwards = TRUE){
  library(tidyverse)
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
  
  
  # STEP 1: Calculate growth rate -------------------------------------------
  x0 <- x
  x <- x %>% 
    arrange(num.from.base) %>% 
    mutate(num.from.tip = rev(num.from.base))
  
  if(method == "auto"){
    # Non-glycine whiskers ----
    if(x$glycine[1] == FALSE){
      ## Assign grouping for the regrowth section 
      # (identify section of the whisker corresponding to the deployment period since all whiskers have a pre-deployment section from the tip)
      # get row number corresponding to start.date by accounting for root length from the tip. 
      l <- x %>% 
        arrange(num.from.tip) %>% 
        mutate(len.from.tip = cumsum(sample.len)) %>%
        pull(len.from.tip)
      rn <- which.min(abs(x$root.length[1] * ratio - l)) + 1 # 0.8 is the ratio of pre-deployment root length to the post-deployment root length. 
      
      x1 <- x %>% 
        arrange(num.from.tip) %>% 
        mutate(rg = 0)
      x1$rg[rn:nrow(x1)] <- seq(1,nrow(x1)-rn+1,1)
      
      ## Calculate growth.rate and resolution ---
      if(x$plucked[1] == TRUE){
        # Plucked non-glycine whiskers
        x <- x1 %>% mutate(res = round(difftime(end.date,start.date, units = 'days')/(original.length - root.length * ratio), 2),
                           growth.rate = (original.length - root.length)/as.integer(difftime(end.date, start.date, units = 'days')))
      } else {
        # Cut non-glycine whiskers
        x <- x1 %>% mutate(res = round(difftime(end.date,start.date, units = 'days')/original.length, 2),
                           growth.rate = original.length/as.integer(difftime(end.date, start.date, units = 'days')))
      }
    } else { 
      ## Glycine seals, plucked ----
      # filter out only regrowth portion based on glycine values. 
      x <- x %>% 
        arrange(num.from.tip) %>% 
        mutate(rg = cumsum(has.glycine))
      
      x1 <- x %>% filter(rg > 0)
      
      # calculate growth.rate and temporal resolution
      x1 <- x1 %>% mutate(res = round(difftime(end.date, start.date, units = 'days')/sum(sample.len), 2), 
                          growth.rate = sum(sample.len)/as.integer(difftime(end.date, start.date, units = 'days')))
      x <- x %>% mutate(res = x1$res[1], growth.rate = x1$growth.rate[1])
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
  


  # STEP 2: Estimate dates for each whisker segment -------------------------
  ### ------ OPTION 1: Starting from end.date (base) ------ ###
  if(estimateBackwards){
    x <- x %>%
      arrange(num.from.base) %>%
      # calculate cumulative length from base
      mutate(len.from.base = cumsum(sample.len)) %>%
      # estimate dates started from the end.date (newest info)
      arrange(len.from.base) %>%
      mutate(date.est = end.date - (lag(len.from.base) * res))
    x$date.est[1] <- x$end.date[1]
  } else if(!estimateBackwards){
    
    ### ------ OPTION 2: Starting from (deployment) start.date ------ ###
    ## First remove predeployment segments and estimate deployment segments
    x1 <- x %>% 
      filter(rg > 0) %>% 
      arrange(num.from.tip) %>% 
      mutate(len.from.start = cumsum(sample.len)) %>%
      mutate(date.est = start.date + (lag(len.from.start) * res))
    
    ## THen estimates dates for segments before start.date (i.e. predeployment)
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
        arrange(num.from.tip)
      
      x2$len.from.start[is.na(x2$date.est)] <- x2$sample.len[is.na(x2$date.est)]
      x2$date.est[is.na(x2$date.est)] <- x2$start.date[1]
      
      ## Combine deployment segments and predeployment segments 
      x <- bind_rows(x2, x1) %>% filter(!is.na(date.est))
    } else {
      x1$date.est[1] <- x1$start.date[1]
      x <- x1
    }
  }
  
  return(x %>% arrange(num.from.base))
}