#=========================================================================================================
#Structure the data
#=========================================================================================================

# Select focal months of study (early sampling was only in these months)
data <- H_killi[which(H_killi$sampling.month%in%c(4,5,6)),]

data$stream <- as.factor(data$stream)
data$sex <- as.factor(data$sex)
data$Habitat <- as.factor(data$Habitat)
data <- data[-which(data$individual_id%in%unique(data[which(data$Dead==1),'individual_id'])),]

data$count <- 1
data <- data[which(data$Habitat=='C' | data$Habitat=='I'),]

data$duplicate<-duplicated(data[,c('stream','individual_id','sampling'),])

data <- data[which(data$duplicate==FALSE),]

data$habitat2 <- as.character(data$Habitat)
data$habitat2 <- as.character(data$Habitat)

data$before_after <- NA
data[which(data$sampling >= 86 & data$sampling<=124),'before_after'] <- '1'
data[which(data$sampling>=170 ),'before_after'] <- '2'



surv <- reshape2::dcast(data,formula=stream  + habitat2 + individual_id + before_after  ~ sampling,value.var='count',fill=0)

data <- data[order(data$sampling),]



#Check for duplicates
surv[grep("2", surv[,c(5:(dim(surv)[2]))]),]


streams <- c("CA", "TY")

survival <- surv[1:4]
survival$captures <-NA 

# create the capture histories
for(i in 1:length(surv$individual_id)){
  
  survival$captures[i] <- paste(surv[i, 5:25], collapse="")
}


survival <- as.data.table(survival)
names(survival) <- c('stream','habitat','individual_id','before_after','captures')

#This removes any individuals that were never caught during these periods

never <- paste(rep("0", nchar(survival$captures[1])), collapse = "")
survival <- survival[captures!=never]

survival$stream <- as.factor(survival$stream)
survival$habitat <- as.factor(survival$habitat)
survival$before_after <- as.factor(survival$before_after)

names(survival) <- c('stream','habitat', 'individual_id', 'before_after','ch')
survival <- survival[,c('ch','stream','habitat', 'before_after')]




#=========================================================================================================
# Run the analyses and perform model averaging
#=========================================================================================================

  for (s in streams){
    
    surv.to.run <- survival[stream==s,]
    
    cap.dates <- ddply(data[which(data$stream==s),],c('stream','sampling'),summarise,cap.date=mean(capture_date,na.rm=TRUE))
    cap.dates$int <- c(cap.dates$cap.date[2:length(cap.dates$cap.date)],NA) - cap.dates$cap.date  
    
    cap.ints <- as.numeric(cap.dates[which(cap.dates$stream==s),'int'])/30
    
    cap.ints <- round(cap.ints,digits = 2)
    cap.ints<- cap.ints[1:length(cap.ints)-1]
      
    dp <- process.data(surv.to.run,model="POPAN",groups=c('habitat'),time.intervals=cap.ints)
    ddl <- make.design.data(dp)
    
 #  Need to specify time periods in the design
    #lets take a look at the time periods
    ddl$Phi$time
    
    # make new design variables for the submodels Phi, p, pent.
    ddl$Phi$period <- NA
    ddl$p$period <- NA
    ddl$pent$period <- NA

    # Define the time periods
    
    # 2015-2016
    ddl$Phi$period[ddl$Phi$Time<20] <- "p1"
    ddl$p$period[ddl$p$Time<20]  <- "p1"
    ddl$pent$period[ddl$pent$Time<20]  <- "p1"
    # 2017-2018
    ddl$Phi$period[ddl$Phi$Time >20 & ddl$Phi$Time <50  ] <- "p2"
    ddl$p$period[ddl$p$Time >20 & ddl$p$Time <50]  <- "p2"
    ddl$pent$period[ddl$pent$Time>20 & ddl$pent$Time<50]  <- "p2"
    #2022-2024
    ddl$Phi$period[ddl$Phi$Time>50] <- "p3"
    ddl$p$period[ddl$p$Time>50]  <- "p3"
    ddl$pent$period[ddl$pent$Time>50]  <- "p3"
    
    ddl$Phi$period <- as.factor(ddl$Phi$period)
    ddl$p$period <- as.factor(ddl$p$period)
    ddl$pent$period <- as.factor(ddl$pent$period)
     
    # models. 
    
  # Phi: survival
    Phi.hab <- list(formula=~habitat)
    Phi.period <- list(formula=~period)
    Phi.hab.period <- list(formula=~habitat+period)
    Phi.habxperiod <- list(formula=~habitat*period)
     
   # p: probability of capture
    p.hab <- list(formula=~habitat)
    p.period <- list(formula=~period)
    p.habperiod <- list(formula=~habitat+period)
    p.habxperiod <- list(formula=~habitat*period)
    
    # pent: probability of entry into population
     pent.hab <- list(formula=~habitat)
     pent.period <- list(formula=~period)
     pent.habperiod <- list(formula=~habitat+period)
     pent.habxperiod <- list(formula=~habitat*period)
     
     # N: super-population size
     N.hab <- list(formula=~habitat)

      cml <- create.model.list("POPAN")
      
      # this is the line of code that does all the heavy lifting and actually fits all the models
      AIC.results <- mark.wrapper.parallel(cml,data=dp,ddl=ddl,delete=TRUE, cpus = 3)

      # model averaged values of Phi

      Phi.average <- model.average(AIC.results, "Phi", alpha=0.025, vcv=TRUE )
      Phi.mean <- as.data.table(Phi.average$estimates)
      
      
      if(s=="CA"){
        Phi.mean <- Phi.mean[(time==1 | time==25.94 | time==86.74) ]
        
         Phi.mean$period <- ifelse(Phi.mean$time==1, 1, 
                                  ifelse(Phi.mean$time==25.94, 2, 3))
     
         }
      
      if(s=="TY"){
        Phi.mean <- Phi.mean[time==1 | time==25.80 | time==86.63] 
        Phi.mean$period <- ifelse(Phi.mean$time==1, 1, 
                                  ifelse(Phi.mean$time==25.80, 2, 3))
      }  
      
      Phi.mean <- Phi.mean[1:6]
      Phi.mean$par <- "Phi"
      Phi.mean$stream <- s
      Phi.vcv <- Phi.average$vcv.real[Phi.mean$par.index, Phi.mean$par.index]

      Phi.average$estimates <- Phi.mean$estimate
      Phi.average$se <- Phi.mean$se
      Phi.average$vcv <- Phi.vcv
      Phi.average$vcv.real <-NA

      results.Phi <- data.table(stream = Phi.mean$stream,
                            habitat = Phi.mean$habitat,
                            period = Phi.mean$period,
                            parameter = Phi.mean$par,
                            estimate = Phi.mean$estimate,
                            se = Phi.mean$se,
                            lcl = Phi.mean$lcl,
                            ucl = Phi.mean$ucl)


      # Extract model averaged derived values (monthly pop size N and births B)

      # weighting of models by AIC score 
      weight <- AIC.results$model.table$weight 

       # number of models
      n.models <- length(weight)

        # number of estimates for population size
      n.estimates <-6

      # matrix to store estimates of pop size
      N.estimate <- matrix(NA,ncol=n.estimates,
                           nrow=n.models)
      # list of matrices to store vcv
      N.vcv <- lapply(X  = 1:n.models, FUN = matrix, data=NA,ncol=n.estimates,
                     nrow=n.estimates)

      # number of estimates of births
      b.estimates <-6
      # matrix to store estimates of births
      B.estimate <- matrix(NA,ncol=b.estimates,
                                nrow=n.models)
      # list of matrices to store vcv
      B.vcv <- lapply(X  = 1:n.models, FUN = matrix, data=NA,ncol=b.estimates,
                      nrow=b.estimates)


      for(i in 1:n.models)
      {
        model.numbers <- as.numeric(row.names(AIC.results$model.table))
        N <- AIC.results[[model.numbers[i]]]$results$derived$`N Population Size`$estimate
       # average over each month
        N.est <- c(mean(N[1:6]), mean(N[7:12]), mean(N[13:21]), 
                   mean(N[22:27]), mean(N[28:33]), mean(N[34:42]))
        
       # add estimates to table
         N.estimate[i,] <- N.est
       
        # var-cov matrix
        Nv <- AIC.results[[model.numbers[i]]]$results$derived.vcv$`N Population Size`
        # average over each month
        
        Nvc <- c(mean(Nv[1:6, 1:6]), mean(Nv[1:6, 7:12]), mean(Nv[1:6, 13:21]), mean(Nv[1:6, 22:27]), mean(Nv[1:6, 28:33]), mean(Nv[1:6, 34:42]),
                 mean(Nv[7:12, 1:6]), mean(Nv[7:12, 7:12]), mean(Nv[7:12, 13:21]), mean(Nv[7:12, 22:27]), mean(Nv[7:12, 28:33]), mean(Nv[7:12, 34:42]),
                 mean(Nv[13:21, 1:6]), mean(Nv[13:21, 7:12]), mean(Nv[13:21, 13:21]), mean(Nv[13:21, 22:27]), mean(Nv[13:21, 28:33]), mean(Nv[13:21, 34:42]),
                 mean(Nv[22:27, 1:6]), mean(Nv[22:27, 7:12]), mean(Nv[22:27, 13:21]), mean(Nv[22:27, 22:27]), mean(Nv[22:27, 28:33]), mean(Nv[22:27, 34:42]),
                 mean(Nv[28:33, 1:6]), mean(Nv[28:33, 7:12]), mean(Nv[28:33, 13:21]), mean(Nv[28:33, 22:27]), mean(Nv[28:33, 28:33]), mean(Nv[28:33, 34:42]),
                 mean(Nv[34:42, 1:6]), mean(Nv[34:42, 7:12]), mean(Nv[34:42, 13:21]), mean(Nv[34:42, 22:27]), mean(Nv[34:42, 28:33]), mean(Nv[34:42, 34:42]))
                 
                 
        # add vcv to list
        N.vcv[[i]]<- matrix(Nvc, nrow=6)

         B <- AIC.results[[model.numbers[i]]]$results$derived$`B Net Births`$estimate
         # average over each between-month period
         B.est <- c(mean(B[1:5]), mean(B[6:11]), mean(B[12:20]), 
                    mean(B[21:25]), mean(B[26:31]), mean(B[32:40]))

         # add estimates to table
         B.estimate[i,] <- B.est #AIC.results[[model.numbers[i]]]$results$derived$`B Net Births`$estimate

         Bv <- AIC.results[[model.numbers[i]]]$results$derived.vcv$`B* Gross Births`

         # average over each between-month period
         
         Bvcv <- c(mean(Bv[1:5, 1:5]), mean(Bv[1:5, 6:11]), mean(Bv[1:5, 12:20]), mean(Bv[1:5, 21:25]), mean(Bv[1:5, 26:31]), mean(Bv[1:5, 32:40]),
                  mean(Bv[6:11, 1:5]), mean(Bv[6:11, 6:11]), mean(Bv[6:11, 12:20]), mean(Bv[6:11, 21:25]), mean(Bv[6:11, 26:31]), mean(Bv[6:11, 32:40]),
                  mean(Bv[12:20, 1:5]), mean(Bv[12:20, 6:11]), mean(Bv[12:20, 12:20]), mean(Bv[12:20, 21:25]), mean(Bv[12:20, 26:31]), mean(Bv[12:20, 32:40]),
                  mean(Bv[21:25, 1:5]), mean(Bv[21:25, 6:11]), mean(Bv[21:25, 12:20]), mean(Bv[21:25, 21:25]), mean(Bv[21:25, 26:31]), mean(Bv[21:25, 32:40]),
                  mean(Bv[26:31, 1:5]), mean(Bv[26:31, 6:11]), mean(Bv[26:31, 12:20]), mean(Bv[26:31, 21:25]), mean(Bv[26:31, 26:31]), mean(Bv[26:31, 32:40]),
                  mean(Bv[32:40, 1:5]), mean(Bv[32:40, 6:11]), mean(Bv[32:40, 12:20]), mean(Bv[32:40, 21:25]), mean(Bv[32:40, 26:31]), mean(Bv[32:40, 32:40]))
         
        # add vcv to list
         B.vcv[[i]] <- matrix(Bvcv, nrow=6)#AIC.results[[model.numbers[i]]]$results$derived.vcv$`B Net Births`

         # tidy up
         rm(model.numbers, N, N.est, Nv, Nvc, B, B.est, Bv, Bvcv)
      }
      # compute model averaged estimates, se, and VCV
      N.average<- model.average(list(estimate=N.estimate,weight=weight,vcv=N.vcv))
      B.average<- model.average(list(estimate=B.estimate,weight=weight,vcv=B.vcv))

      # Now we need to convert N to population density,
      # and B to per capita recruitment

     #  # population density (fish per metre length of stream)
       streamlength<- stream_data[stream==s, stream.length]
       streamlength.mat <- streamlength %*% t(streamlength)

       pop.dens <- results.Phi
       pop.dens$parameter <- "pop.dens"
       pop.dens$estimate <- N.average$estimate / streamlength
       pop.dens$se <- N.average$se / streamlength
       pop.dens$lcl <-  pop.dens$estimate - 1.96*pop.dens$se 
       pop.dens$ucl <- pop.dens$estimate + 1.96*pop.dens$se 
       
 # per capita recruitment
      recruitment <- results.Phi
      recruitment$parameter <- "babies"
      recruitment$estimate <- B.average$estimate / N.average$estimate
      recruitment$se <- B.average$se / N.average$estimate
      recruitment$lcl <- recruitment$estimate - 1.96*recruitment$se 
      recruitment$ucl <- recruitment$estimate + 1.96*recruitment$se
      
    # add them into the main results
      
      results <- rbind(results.Phi, pop.dens, recruitment)

    if(s=='TY'){
        AIC.results.TY <- AIC.results
        # # save AIC scores and weighting of different models
        write.csv(AIC.results.TY$model.table,file='./output/results/AICTY.csv')
        # save all models from model selection process
        save(AIC.results.TY, file="./output/results/POPAN-models-TY.Rda")

    }
    if(s=='CA'){
        AIC.results.CA <- AIC.results
        write.csv(AIC.results.CA$model.table,file='./output/results/AICCA.csv')
        # save all models from model selection process
        save(AIC.results.CA, file="./output/results/POPAN-models-CA.Rda")
    }

      if(s==streams[1]){
        POPAN.parameters.all <- results
      }
     
      if(s!=streams[1]){
        POPAN.parameters.all <- rbind(POPAN.parameters.all,results)
        # save the combined parameters for both streams
        write.csv(POPAN.parameters.all, file="./output/results/POPAN-parameters.csv")

      }
    
    print(s)
    
    # tidy up
     rm(Phi.hab, Phi.period, Phi.hab.period, Phi.habxperiod, 
        p.hab, p.period, p.habperiod, p.habxperiod, 
        pent.hab, pent.period, pent.habperiod, pent.habxperiod, N.hab
        )
    
  }


