# SPAS input file formatting function

SPAS.input.file <- function(appdat, recdat, app.strata, combine.app,
                            rec.strata, combine.rec,
                            sex.correction, name, sex.include = c("M", "F")){ 
  
  # If combine.app is TRUE, this will unite the strata for subsequent grouping
  if(combine.app == TRUE) combine.by.a <- enquo(app.strata)
  if(combine.app == TRUE) appdat <- appdat %>% unite("app.strata", !!combine.by.a)
  if(combine.app == TRUE) app.strata = quo(app.strata)
  # If combine.app is FALSE, this will just wrap the strata in qosure o that subsequent functions can read it as a column name
  if(combine.app == FALSE) app.strata <- enquo(app.strata)
  
  # If combine.rec is TRUE, this will unite the strata for subsequent grouping
  if(combine.rec == TRUE) combine.by.r <- enquo(rec.strata)
  if(combine.rec == TRUE) recdat <- recdat %>% unite("rec.strata", !!combine.by.r)
  if(combine.rec == TRUE) appdat <- appdat %>% unite("rec.strata", !!combine.by.r)
  if(combine.rec == TRUE) rec.strata = quo(rec.strata)
  # If combine.rec is FALSE, this will just wrap the strata in qosure o that subsequent functions can read it as a column name
  if(combine.rec == FALSE) rec.strata <- enquo(rec.strata)
  
  correction = sex.correction
  
  # create sex corrected applications~app strata  
  sex.cor <- filter(appdat, Application.Sex %in% c("M", "F")) %>% 
    group_by(!!app.strata) %>% # creating a row for each app strata
    summarise(applied = n()) %>% # n() plus group by does the same thing as count()
    mutate(Tagged.Males = applied*(correction),
           Tagged.Females = -Tagged.Males,
           Tagged.Jacks = 0) # the sex correction per strata
  
  # This function is two sets of the "same" code, one for Males and one for Females
  # Male code is denoted by 
  # Female code is denoted by f.
  # Only the first set of steps is annotated
  
  for(i in sex.include){
    
    # quosures are used by TidyVerse functions. Only went this route because I had already done everything in tidyverse.
    # the below quosures are used throughout the rest of the loop
    unt.sex <- switch(i, "M" = quo(Untagged.Males), "F" = quo(Untagged.Females), "J" = quo(Untagged.Jacks))
    tag.sex <- switch(i, "M" = quo(Tagged.Males), "F" = quo(Tagged.Females), "J" = quo(Tagged.Jacks))
    SEX <- switch(i, "M" = "Male", "F" = "Female", "J" = "Jack")
    
    # add sex corrected applications~app strata 
    tagged.sex.cor <- pull(sex.cor, !!tag.sex)
    
    cor.app <- filter(appdat, Application.Sex == i) %>% 
      group_by(!!app.strata) %>% # creating a row for each app strata
      summarise(applied = n()) %>% # n() plus group by does the same thing as count()
      mutate(applied = applied + tagged.sex.cor, # adding the sex correction
             pr = applied/sum(applied)) # proportion applied in each strata
    
    # create sum of Untagged fish recovered per recovery strata
    unt.rec <- recdat %>% 
      group_by(!!rec.strata) %>%
      summarise(untagged.rec = sum(!!unt.sex))
    
    # create column for total recoveryies and proportion recovered in each strata that will be used for redistributing secondary marked fish to application strata proportional to the recovery strata they were found in
    total.rec <- recdat %>% 
      group_by(!!rec.strata) %>%
      summarise(total.rec = sum(!!tag.sex + !!unt.sex)) %>%
      mutate(pr = total.rec/sum(total.rec)) # proportion recovered in each strata
    
    # Check to see if there were NO marked recoveries in any of the strata
    CHECK.marked <- merge(total.rec, unt.rec) %>%
      mutate(marked = total.rec - untagged.rec)
    
    if(any(CHECK.marked$marked <= 0)) cat(paste0("No marked ", SEX, " recoveries in strata \n Input file for ", SEX, " not created"))
    # if there are no marked recoveries in one of the recovery strata, the loop stops here and produces the above message
    if(any(CHECK.marked$marked > 0)){ # if all good, the loop continues
      
      # t by s matrix of number of tags recovered in strata t by application strata s
      t.rec <- appdat %>%
        filter(Recovery.Sex == i) %>%
        group_by(!!app.strata, !!rec.strata) %>%
        summarise(tag.recoverd = n()) %>% # this creates a long form summary
        spread(!!rec.strata, tag.recoverd) # this changes it to the wide form t by s matrix
      if(combine.app == TRUE & t.rec[nrow(t.rec),1] == "NA_NA") t.rec[nrow(t.rec),1] <- NA
      
      add <- setdiff(cor.app[,1], data.frame(app.strata = na.omit(t.rec[,1])))
      num <- nrow(add)
      if(num >0) t.rec[(nrow(t.rec)+1):(nrow(t.rec)+num), 1] <- add
      t.rec <- arrange(t.rec, !!app.strata) # only works if its dates or the app strata are alphabetical
      
      t.rec <- as.data.frame(t.rec) # needed to do this because funky things were happening when it was a tibble
      
      # need to fill in if any recovery spot in the app strata had a 0, all of the subsetting is to ignore if there are NAs in the app_strata column as it is used as a check to apply secondary mark only recoveries
      t.rec[ , 2:ncol(t.rec)][is.na(t.rec[ , 2:ncol(t.rec)])] <- 0
      
      
      # if there are secondary only recoveries the next three blocks will run
      
      # takes the proportion column from cor.app and uses matrix multiplication to distribute the number of secondary recoveries from each t strata to a t by s matrix
      if(is.na(t.rec[nrow(t.rec),1])) secondary.add <- as.matrix(cor.app[,3]) %*% 
        as.matrix(t.rec[is.na(t.rec[,1]) , 2:ncol(t.rec)])
      
      # if there are any NA's because there were no secondary only recoveries in a particular strata, this replaces them with a 0 (I don't know why I had to make an equivelant number of 0s to NAs...)  
      if(is.na(t.rec[nrow(t.rec),1])) secondary.add[is.na(secondary.add)] <- 
        rep(0, length(secondary.add[is.na(secondary.add)]))
      
      # adds the secondary correction to the t by s matrix of tag recoveries
      if(is.na(t.rec[nrow(t.rec),1])) t.rec[!is.na(t.rec[,1]), 
                                            2:ncol(t.rec)] <- t.rec[!is.na(t.rec[,1]), 
                                                                    2:ncol(t.rec)] + secondary.add
      
      # assembles new R SPAS input matrix & file
      spasR <- cbind( t.rec[!is.na(t.rec[,1]), -1],  cor.app[ , 2] - rowSums(t.rec[!is.na(t.rec[,1]), -1]))
      spasR <- as.matrix(rbind(spasR, c(t(unt.rec[ , -1]), 0)))
      
      write.table(spasR, paste0("SPAS_R_", name, "_input_", SEX, "s.csv"), 
                  row.names = F, col.names = F, sep = ",")
      write.table(cor.app[,1], paste0("SPAS_R_", name, "_for_", SEX, "s_ROW_STRATA_NAMES.csv"), sep = ",")
      
      cat(paste0("\n", "SPAS_R_", name, " input and strata names has been created for ", SEX, "s \n \n"))
      cat(paste0("Sex corrected values for comparison with ", SEX, 
                 " PopEstimatorPRT tab in WB5PopEstimator.xlsx \n * Note that only Applied may be different due to sex correction \n \n"))
      cat(SEX)
      cat("Tags Applied:", sum(spasR[-nrow(spasR), ]), "\n")
      cat("Tagged Recoveries:", sum(spasR[-nrow(spasR), -ncol(spasR)]), "\n")
      cat("Total Recoveries:", sum(spasR[,-ncol(spasR)]), "\n")
      cat("\n", "------------------------------------------------------------", "\n")
    } # end of sex.include loop
  } # end of IF statement if no tags recovered in one or more recovery strata
} # end of SPAS.input.file function


# define reporting and removal functions
get.models <- function(pattern){
  # get the list of models according to the pattern specified
  #browser()
  model.list <- mget( ls(envir=globalenv())[grepl(pattern,ls(envir=globalenv()))], envir=globalenv())
  model.list
}


make.report <- function(model.list){
  # make a little report
  report <- plyr::ldply(model.list, function(x){
    #browser()
    data.frame(#version=x$version,
      date   = as.Date(x$date),
      model.id         = x$model.info$model.id,
      AICc             = round(x$model.info$AICc, 3),
      Nhat             = round(x$est$real$N),
      Nhat.se          = round(x$se $real$N),
      s.a.pool         =-1+nrow(x$fit.setup$pooldata),
      t.p.pool         =-1+ncol(x$fit.setup$pooldata),
      logL.cond        = x$model.info$logL.cond,
      np               = x$model.info$np,
      chapman.est      = round(x$est$N.Chapman, 0),
      kappa.pooled     = round(x$kappa.after.lp, 0),
      gof.chisq        = round(x$gof$chisq, 1),
      gof.df           = x$gof$chisq.df,
      gof.p            = round(x$gof$chisq.p, 3)
    )
    
  })
  report <- mutate(report, deltaAIC = round(AICc - min(AICc)), 3) %>% arrange(deltaAIC)
  report <- report[, c(1:3,16,4:15)]
  report
}