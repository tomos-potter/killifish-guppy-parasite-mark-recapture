#=========================================================================================================
#Structure the data
#=========================================================================================================

# Select the months of study
data <- H_killi[which(H_killi$sampling.month%in%c(4,5,6)),]

streams <- c("CA", "TY")


data$stream <- as.factor(data$stream)
data$sex <- as.factor(data$sex)
data$Habitat <- as.factor(data$Habitat)
data <- data[-which(data$individual_id%in%unique(data[which(data$Dead==1),'individual_id'])),]

data <- data[which(data$Habitat=='C' | data$Habitat=='I'),]

data$duplicate<-duplicated(data[,c('stream','individual_id','sampling'),])


data$habitat2 <- as.character(data$Habitat)

#Remove duplicated captures
data$duplicate<-duplicated(data[,c('stream','individual_id','sampling'),])
data <- data[-which(data$duplicate==TRUE),]

# define size classes
data[which(data$TL>15 & data$TL<=35),'size'] <- 1
data[which(data$TL>35 & data$TL<=55),'size'] <- 2
data[which(data$TL>55 ),'size'] <- 3

#=========================================================================================================
#Create capture histories
#=========================================================================================================

data$habitat2 <- as.character(data$Habitat)

data$before_after <- NA
data[which(data$sampling >= 86 & data$sampling<=124),'before_after'] <- 'B'
data[which(data$sampling>=170 ),'before_after'] <- 'A'

surv <- reshape2::dcast(data,formula=stream  + habitat2 + before_after + individual_id ~ sampling,value.var='size',fill=0)
data <- data[order(data$sampling),]


#Check for duplicates
surv[grep("2", surv[,c(5:(dim(surv)[2]))]),]


streams <- c("CA", "TY")

survival <- surv[1:4]
survival$captures <-NA 

for(i in 1:length(surv$individual_id)){
  
  survival$captures[i] <- paste(surv[i, 5:25], collapse="")
}


survival <- as.data.table(survival)
names(survival) <- c('stream','habitat','before_after','individual_id','captures')

#This removes any individuals that were never caught during these periods

never <- paste(rep("0", nchar(survival$captures[1])), collapse = "")
survival <- survival[captures!=never]

survival$stream <- as.factor(survival$stream)
survival$habitat <- as.factor(survival$habitat)
survival$before_after <- as.factor(survival$before_after)

names(survival) <- c('stream','habitat','before_after', 'individual_id','ch')
survival <- survival[,c('ch','stream','habitat', 'before_after')]


#=========================================================================================================
#Run Analyses
#=========================================================================================================

for (s in streams){
  
  surv.to.run <- survival[stream==s,]
  
  # define the capture intervals
  cap.dates <- ddply(data[which(data$stream==s),],c('stream','sampling'),summarise,cap.date=mean(capture_date,na.rm=TRUE))
  cap.dates$int <- c(cap.dates$cap.date[2:length(cap.dates$cap.date)],NA) - cap.dates$cap.date  
  
  cap.ints <- as.numeric(cap.dates[which(cap.dates$stream==s),'int'])/30
  
  cap.ints <- round(cap.ints,digits = 2)
  cap.ints<- cap.ints[1:length(cap.ints)-1]
  
  # initiate the models
  dp <- process.data(surv.to.run,model="Multistrata",groups=c('habitat'), time.intervals=cap.ints)
  ddl <- make.design.data(dp)
  
  #  Need to specify time periods in the design
  #lets take a look at the time periods
  ddl$S$time
  
  # make new design variables for the submodels Phi, p, pent.
  ddl$S$period <- NA
  ddl$p$period <-NA
  ddl$Psi$period <-NA

  # Define the time periods
  # 2015-2016
  ddl$S$period[ddl$S$Time<20] <- "p1"
  ddl$p$period[ddl$p$Time<20]  <- "p1"
  ddl$Psi$period[ddl$Psi$Time<20]  <- "p1"
  # 2017-2018
  ddl$S$period[ddl$S$Time >20 & ddl$S$Time <50  ] <- "p2"
  ddl$p$period[ddl$p$Time >20 & ddl$p$Time <50]  <- "p2"
  ddl$Psi$period[ddl$Psi$Time>20 & ddl$Psi$Time<50]  <- "p2"
  # 2022-2024
  ddl$S$period[ddl$S$Time>50] <- "p3"
  ddl$p$period[ddl$p$Time>50]  <- "p3"
  ddl$Psi$period[ddl$Psi$Time>50]  <- "p3"
  
  
  # fix some transition probabilities
  # we are only modelling transitions from 1 to 2 and from 2 to 3
  ddl$Psi$fix <- NA
  ddl$Psi$fix[ddl$Psi$stratum==1 & ddl$Psi$tostratum==3]<-0
  ddl$Psi$fix[ddl$Psi$stratum==3 & ddl$Psi$tostratum==1]<-0
  ddl$Psi$fix[ddl$Psi$stratum==3 & ddl$Psi$tostratum==2]<-0
  ddl$Psi$fix[ddl$Psi$stratum==2 & ddl$Psi$tostratum==1]<-0

  
  # Models. 

  # survival  
  S.periodXhab <- list(formula=~ habitat*period)
  S.sizeplusperiodXhab <- list(formula=~ habitat*period + stratum)
  S.sizeXhabXperiod <-list(formula=~ (habitat*period) * stratum)

  # recapture probability
  p.habper  <- list(formula=~ habitat*period)
  p.sizePLUShabper <- list(formula=~ habitat*period + stratum)
  p.sizeXhabper <- list(formula=~ (habitat*period) * stratum)

  # probability of transition between sizes
 # Psi.size <- list(formula=~-1 + stratum:tostratum)
  Psi.sizePLUShabper <- list(formula=~ (stratum:tostratum)+(habitat*period), remove.intercept=TRUE)
  Psi.sizeXhabperiod <- list(formula=~ (stratum:tostratum):(habitat*period), remove.intercept=TRUE)

  
  # Runs model selection and spits out results
    cml <- create.model.list("Multistrata")
    AIC.results <- mark.wrapper.parallel(cml,data=dp,ddl=ddl,delete=TRUE,cpus=3)
  # model averaging for S, weighted by AIC
     S.average <- model.average(AIC.results, "S", alpha=0.025, vcv=TRUE, drop=TRUE )
     surv.mean <- as.data.table(S.average$estimates)
     surv.mean <- surv.mean <- surv.mean[(occ==1 | occ==7 |occ==13) & cohort==1, 
                                         c(1:5,12,14,18)]

     surv.mean$period <- ifelse(surv.mean$occ==1, 1, ifelse(surv.mean$occ==7, 2, 3))
     surv.mean$par <- "S"
     surv.mean$stream <- s

   # variance-covariance matrix for model averaged results                            
     surv.vcv <- S.average$vcv.real[surv.mean$par.index, surv.mean$par.index]
     S.average$vcv <- surv.vcv
     S.average$estimates <- surv.mean$estimate
     S.average$vcv.real <-NULL
    
  # record the results in a table 
    results.S <- data.table(stream = surv.mean$stream,
                            habitat = surv.mean$habitat,
                            period = surv.mean$period,
                            parameter = surv.mean$par,
                            stratum = surv.mean$stratum,
                            estimate = surv.mean$estimate,
                            se = surv.mean$se,
                            lcl = surv.mean$lcl,
                            ucl = surv.mean$ucl)
    
  # model averaging for Psi, weighted by AIC
     Psi.average <- model.average(AIC.results, "Psi", alpha=0.025, vcv=TRUE  )
     Psi.mean <- as.data.table(Psi.average$estimates)
     Psi.mean <- Psi.mean[(occ==1 | occ==7 |occ==13) & fixed!="Fixed" & occ.cohort==1,
                          c(1:5,12,14:15,19)]
     
     Psi.mean$stratum <- paste(Psi.mean$stratum, Psi.mean$tostratum, sep="to")
     Psi.mean$period <- ifelse(Psi.mean$occ==1, 1, ifelse(Psi.mean$occ==7, 2, 3))
     
     Psi.mean$par <- "Psi"
     Psi.mean$stream <- s
     Psi.mean$par.index <- Psi.mean$par.index - dim(Psi.average$estimates)[1]
     
 # variance-covariance matrix for model averaged results                            
     Psi.average$vcv <- Psi.average$vcv.real[Psi.mean$par.index, Psi.mean$par.index]
     Psi.average$vcv.real <-NULL
     Psi.average$estimates <- Psi.mean$estimate
    
 
     results.Psi <- data.table(stream = Psi.mean$stream,
                             habitat = Psi.mean$habitat,
                             period = Psi.mean$period,
                             parameter = Psi.mean$par,
                             stratum = Psi.mean$stratum,
                             estimate = Psi.mean$estimate,
                             se = Psi.mean$se,
                             lcl = Psi.mean$lcl,
                             ucl = Psi.mean$ucl)
# combine S and Psi results   
     results <- rbind(results.S, results.Psi)
    
     Multistrata.results <- list(results = results,
                                 S = S.average,
                                 Psi = Psi.average)

  if(s=='TY'){

      AIC.results.TY.size <- AIC.results
      # save all models from model selection process
      save(AIC.results.TY.size, file="./output/results/Multistrata-models-TY.Rda")
      # write a table of AIC scores for each model
      write.csv(AIC.results.TY.size$model.table,file='./output/results/AICTYsize.csv')

  }
  if(s=='CA'){
    AIC.results.CA.size <- AIC.results
    save(AIC.results.CA.size, file="./output/results/Multistrata-models-CA.Rda")
    write.csv(AIC.results.CA.size$ model.table,file='./output/results/AICCAsize.csv')
  }
  
     
     if(s==streams[1]){
       Multistrata.parameters.all <- results
     }
     
     if(s!=streams[1]){
       Multistrata.parameters.all <- rbind(Multistrata.parameters.all,results)
       # save the combined parameters for both streams
       write.csv(Multistrata.parameters.all, file="./output/results/Multistrata-parameters.csv")
       
    }

}



