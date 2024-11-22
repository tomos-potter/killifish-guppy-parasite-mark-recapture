#=========================================================================================================
#Structure the data
#=========================================================================================================

# Select focal months of study (early sampling was only in these months)
data <- H_guppy[which(H_guppy$sampling.month%in%c(4,5,6)),]
data <- as.data.table(data)
data$stream <- as.factor(data$stream)
data$sex <- as.factor(data$sex)
data$Habitat <- as.factor(data$Habitat)
#Removes individuals that died during capture
gonners <- data[Dead==1, individual_id]

data <- data[ ! data$individual_id %in% gonners, ]
data$count <- 1
data$habitat2 <- as.character(data$Habitat)
surv <- reshape2::dcast(data,formula=stream  + individual_id ~ sampling,fill=0)

#Removes duplicates
 grep(2,surv[,c(3,14)])
 surv[grep(2,surv$'172'),'172'] <- 1


data$count <- 1
surv <- reshape2::dcast(data,formula=stream  + individual_id ~ sampling,value.var='count',fill=0)

data <- data[order(data$sampling),]

#check for duplicates
surv[grep("2", surv[,c(5:(dim(surv)[2]))]),]
surv[grep("3", surv[,c(5:(dim(surv)[2]))]),]
surv[grep("4", surv[,c(5:(dim(surv)[2]))]),]

streams <- c("CA", "TY")

survival <- surv[1:2]
survival$captures <-NA 

# create the capture histories
for(i in 1:length(surv$individual_id)){
  
  survival$captures[i] <- paste(surv[i, 3:23], collapse="")
}


survival <- as.data.table(survival)

names(survival) <- c('stream','individual_id', 'ch')
survival <- survival[,c('ch','stream')]

# remove fish "caught" twice in one month
survival <- survival[-grep("2", survival$ch),]

#=========================================================================================================
# Run the analyses and perform model averaging
#=========================================================================================================

for (s in streams){
 
  surv.to.run <- survival[which(survival$stream==s ),]
  
  cap.dates <- ddply(data[which(data$stream==s),],c('stream','sampling'),summarise,cap.date=mean(capture_date,na.rm=TRUE))
  cap.dates$int <- c(cap.dates$cap.date[2:length(cap.dates$cap.date)],NA) - cap.dates$cap.date  
  
  cap.ints <- as.numeric(cap.dates[which(cap.dates$stream==s),'int'])/30
  
  cap.ints <- round(cap.ints,digits = 2)
  cap.ints<- cap.ints[1:length(cap.ints)-1]
  
  dp <- process.data(surv.to.run,model="POPAN",groups=c('stream'),time.intervals=cap.ints)
  ddl <- make.design.data(dp)
  
  #  Need to specify time periods in the design
  #lets take a look at the time periods
  ddl$Phi$time
  
  # make new design variables for the submodels Phi, p, pent.
  ddl$Phi$period <- NA
  ddl$p$period <-NA
  ddl$pent$period <-NA

  # Define the time periods
    
    ddl$Phi$period[ddl$Phi$Time<20] <- "p1"
    ddl$p$period[ddl$p$Time<20]  <- "p1"
    ddl$pent$period[ddl$pent$Time<20]  <- "p1"
    
    ddl$Phi$period[ddl$Phi$Time >20 & ddl$Phi$Time <50  ] <- "p2"
    ddl$p$period[ddl$p$Time >20 & ddl$p$Time <50]  <- "p2"
    ddl$pent$period[ddl$pent$Time>20 & ddl$pent$Time<50]  <- "p2"
    
    ddl$Phi$period[ddl$Phi$Time>50] <- "p3"
    ddl$p$period[ddl$p$Time>50]  <- "p3"
    ddl$pent$period[ddl$pent$Time>50]  <- "p3"

  
  # models. 
  Phi.dot  <-  list(formula=~1)
  Phi.period <- list(formula=~period)
  
  p.dot  <-  list(formula=~1)
  p.period <- list(formula=~period)
  
  pent.dot <- list(formula=~1)
  pent.period <- list(formula=~period)
  
  N.dot <- list(formula=~1)

    cml <- create.model.list("POPAN")
    AIC.results <- mark.wrapper.parallel(cml,data=dp,ddl=ddl,delete=TRUE, cpus=3)
   

    # model averaged values of Phi
    Phi.average <- model.average(AIC.results, "Phi", alpha=0.025, vcv=TRUE, drop=TRUE )
    Phi.mean <- as.data.table(Phi.average$estimates)
    Phi.mean <- Phi.mean[par.index==1 | par.index==7 |par.index==13,]
    Phi.mean$par <- "Phi"
    Phi.mean$period <- ifelse(Phi.mean$par.index==1, 1,
                              ifelse(Phi.mean$par.index==7, 2,3))
    Phi.mean$stream <- s
    Phi.vcv <- Phi.average$vcv.real[Phi.mean$par.index, Phi.mean$par.index]

    
    Phi.average$estimates <- Phi.mean$estimate
    Phi.average$se <- Phi.mean$se
    Phi.average$vcv <- Phi.vcv
    
    results.Phi <- data.table(stream = Phi.mean$stream,
                              period = Phi.mean$period,
                              parameter = Phi.mean$par,
                              estimate = Phi.mean$estimate,
                              se = Phi.mean$se,
                              lcl = Phi.mean$lcl,
                              ucl = Phi.mean$ucl)
    
    # Extract model averaged derived values (monthly pop size N and births B)
    
    # weighting of models by AIC score (extract only those with weight>0)
    weight <- AIC.results$model.table$weight[which(AIC.results$model.table$weight>0)]
    
    # number of models
    n.models <- length(weight)
    
    # number of estimates for population size
    n.estimates <-3
    
    # matrix to store estimates of pop size
    N.estimate <- matrix(NA,ncol=n.estimates,
                         nrow=n.models)
    # list of matrices to store vcv
    N.vcv <- lapply(X  = 1:n.models, FUN = matrix, data=NA,ncol=n.estimates,
                    nrow=n.estimates)
    
    # number of estimates of births
    b.estimates <-3
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
      N.est <- c(mean(N[1:6]), mean(N[7:12]), mean(N[13:21]))
      # add estimates to table
      N.estimate[i,] <- N.est
      
      Nv <- AIC.results[[model.numbers[i]]]$results$derived.vcv$`N Population Size`
      # average over each month
      Nvc <- c(mean(Nv[1:6, 1:6]), mean(Nv[1:6, 7:12]), mean(Nv[1:6, 13:21]),
               mean(Nv[7:12, 1:6]), mean(Nv[7:12, 7:12]), mean(Nv[7:12, 13:21]),
               mean(Nv[13:21, 1:6]), mean(Nv[13:21, 7:12]), mean(Nv[13:21, 13:21]))
      # add vcv to list
      N.vcv[[i]]<- matrix(Nvc, nrow=3)
      
      B <- AIC.results[[model.numbers[i]]]$results$derived$`B Net Births`$estimate
      # average over each between-month period
      B.est <- c(mean(B[1:6]), mean(B[7:12]), mean(B[13:20]))
      
      # add estimates to table
      B.estimate[i,] <- B.est #AIC.results[[model.numbers[i]]]$results$derived$`B Net Births`$estimate
      
      Bv <- AIC.results[[model.numbers[i]]]$results$derived.vcv$`B Net Births`
      
      # average over each between-month period
      Bvcv <- c(mean(Bv[1:6, 1:6]), mean(Bv[1:6, 7:12]), mean(Bv[1:6, 13:20]),
                mean(Bv[7:12, 1:6]), mean(Bv[7:12, 7:12]), mean(Bv[7:12, 13:20]),
                mean(Bv[13:20, 1:6]), mean(Bv[13:20, 7:12]), mean(Bv[13:20, 13:20]))

      
      # add vcv to list 
      B.vcv[[i]] <- matrix(Bvcv, nrow=3)#AIC.results[[model.numbers[i]]]$results$derived.vcv$`B Net Births`
      
    }
    # compute model averaged estimates, se, and VCV
    N.average<- model.average(list(estimate=N.estimate,weight=weight,vcv=N.vcv), drop=TRUE)
    B.average<- model.average(list(estimate=B.estimate,weight=weight,vcv=B.vcv), drop=TRUE)
    
    # Now we need to convert N to population density,
    # and B to per capita recruitment
    
    # population density (fish per metre length of stream)
    streamlength<- stream_data[stream==s & period==3 & habitat=="I", stream.length]
    streamlength2 <- streamlength^2
    streamlength.mat <- matrix(rep(streamlength, 9), nrow=3)
    
    pop.dens <- N.average
    pop.dens$estimate <- N.average$estimate / streamlength
    pop.dens$se <- N.average$se / streamlength
    pop.dens$vcv <- N.average$vcv/streamlength.mat
    
    # per capita recruitment
    
    recruitment <- B.average
    recruitment$estimate <- B.average$estimate / N.average$estimate
    recruitment$se <- B.average$se / N.average$estimate
    
    recruitment$vcv <- B.average$vcv / (N.average$estimate %*% t(N.average$estimate))
    
    # add them into the main results
    results.popdens <- results.Phi
    results.popdens$parameter <-"pop.dens"
    results.popdens$estimate <- pop.dens$estimate
    results.popdens$se <- pop.dens$se
    results.popdens$lcl <- pop.dens$estimate - (1.96 * pop.dens$se)
    results.popdens$ucl <- pop.dens$estimate + (1.96 * pop.dens$se)
    
    results.rec <- results.Phi
    results.rec$parameter <-"per.capita.recruitment"
    results.rec$estimate <- recruitment$estimate
    results.rec$se <- recruitment$se
    results.rec$lcl <- recruitment$estimate - (1.96 * recruitment$se)
    results.rec$ucl <- recruitment$estimate + (1.96 * recruitment$se)
    
    results <- rbind(results.Phi, results.popdens, results.rec)
    
    
    POPAN.results <- list(results = results, 
                          Phi = Phi.average, 
                          N = N.average, 
                          B = B.average, 
                          pop.dens = pop.dens, 
                          recruitment = recruitment)
    
    
    if(s=='TY'){
      
      
      AIC.results.TY <- AIC.results
      # save AIC scores and weighting of different models
      write.csv(AIC.results.TY$ model.table,file='./output/results/G-AICTY.csv')
      # save the model results list
      save(AIC.results.TY, file = './output/results/G-POPAN-models-TY.Rda')
      
    }
    if(s=='CA'){
      
      
      AIC.results.CA <- AIC.results
      write.csv(AIC.results.CA$ model.table,file='./output/results/G-AICCA.csv')
      # save the model results list
      save(AIC.results.CA, file = './output/results/G-POPAN-models-CA.Rda')
      
    }
    


  if(s==streams[1]){
    surv.results.all <- results
  }
  
  if(s!=streams[1]){
    surv.results.all <- rbind(surv.results.all,results)
    
    write.csv(surv.results.all, file="./output/results/G-POPAN-parameters.csv")
  }
  
 
  
  print(s)
}

