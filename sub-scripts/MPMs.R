# Stage-structured matrix projection models
#=========================================================================================================
# Load the parameters and format
#=========================================================================================================

pars <- fread("./output/results/Multistrata-parameters.csv")
pars2 <- fread("./output/results/POPAN-parameters.csv")

pars$par.name <- paste(pars$parameter, pars$stratum, sep=".")
pars$treatment <- paste(pars$stream, pars$habitat, pars$period, sep="-")
pars2$treatment <- paste(pars2$stream, pars2$habitat, pars2$period, sep="-")

wide_par <- pivot_wider(pars[,c("stream","habitat","period","par.name", "treatment", "estimate", "se")], 
                   names_from = par.name,
                   values_from = c("estimate", "se"))

wide_par2 <- pivot_wider(pars2[,c("stream","habitat","period","parameter","treatment", "estimate", "se")], 
                         names_from = parameter,
                         values_from = c("estimate", "se"))

parameters <- merge(x=wide_par, y=wide_par2, by=c("treatment", "stream",
                                                  "habitat", "period"))

#=========================================================================================================
# FUNCTIONS
#=========================================================================================================

# function that constructs the matrix projection model

three_stage_matrix <- function(surv1, surv2, surv3,growth12, growth23,fec1, fec2, fec3) {
  
  stay1 = surv1*(1-growth12)
  grow12 = surv1*growth12
  stay2 = surv2*(1-growth23)
  grow23 = surv2*(growth23)
  stay3 = surv3
 
  rec2 = fec2 * surv1
  rec3 = fec3 * surv1
  
  mat <- matrix(data=c(stay1,     rec2,     rec3,
                       grow12,  stay2,    0,    
                       0,         grow23, stay3), ncol=3, byrow=TRUE)
  return(mat)
}

# function that propagates parameter uncertainty by
# using values drawn from the mean and se of estimates,
# then calculates rv1 for each,
# returns rv1 mean, se, and 95% confidence intervals

mean.model <- function(p, it){
  
  set.seed(1987)
  
  # to store results
  results <- data.table(stream = p["stream"],
                       hab = p["habitat"],
                       period = p["period"],
                       it = it)
  
  # to store estimates of rv1
  squid <- data.table(stream = p["stream"],
                      hab = p["habitat"],
                      period = p["period"],
                      it = 1:it,
                      surv1 = NA,
                      surv2 = NA,
                      surv3 = NA,
                      growth12 = NA,
                      growth23 = NA,
                      fec2 = NA,
                      lambda = NA, # population growth rate at equilibrium
                      ss1 = NA,   # stable stage structure: proportion of size class 1
                      ss2 = NA,   # stable stage structure: proportion of size class 2
                      ss3 = NA,   # stable stage structure: proportion of size class 3
                      rv1 = NA,   # reproductive value size class 1
                      rv2 = NA,   # reproductive value size class 2
                      rv3 = NA)   # reproductive value size class 3
  for (i in 1:it){
    
    #generate random draws from the probability distribution of each parameter
    squid$surv1[i] = betaval(p$estimate_S.1, sd = p$se_S.1)
    squid$surv2[i] = betaval(p$estimate_S.2, sd = p$se_S.2)
    squid$surv3[i] = betaval(p$estimate_S.3, sd = p$se_S.3)
    squid$growth12[i] = betaval(p$estimate_Psi.1to2, sd = p$se_Psi.1to2)
    squid$growth23[i] = betaval(p$estimate_Psi.2to3, sd = p$se_Psi.2to3)
    squid$fec2[i] = lnorms(1, (p$estimate_babies), (p$se_babies)^2)*0.5
    squid$fec3[i] = lnorms(1, (p$estimate_babies), (p$se_babies)^2)*1.5
    
    # build the MPM
    x <- three_stage_matrix(squid$surv1[i], squid$surv2[i], squid$surv3[i],
                            squid$growth12[i], squid$growth23[i], 
                            squid$fec1[i], squid$fec2[i], squid$fec3[i])
  
   # calculate lambda 
   squid$lambda[i] = lambda(x)
   
   #calculate stable stage distribution
   ssd <- stable.stage(x)
   
   squid$ss1[i] = ssd[1]
   squid$ss2[i] = ssd[2]
   squid$ss3[i] = ssd[3]
   
   #calculate reproductive values
   rv <- reproductive.value(x)
   squid$rv1[i] = rv[1]
   squid$rv2[i] = rv[2]
   squid$rv3[i] = rv[3]
    
  }

  results$lambda.mean = mean(squid$lambda)
  results$lambda.se = sd(squid$lambda)
  results$lambda.lcl = results$lambda.mean - 1.96*results$lambda.se
  results$lambda.ucl = results$lambda.mean + 1.96*results$lambda.se
  
  results$ss1.mean = mean(squid$ss1)
  results$ss1.se = sd(squid$ss1)
  results$ss1.lcl = results$ss1.mean - 1.96*results$ss1.se
  results$ss1.ucl = results$ss1.mean + 1.96*results$ss1.se
  
  results$ss2.mean = mean(squid$ss2)
  results$ss2.se = sd(squid$ss2)
  results$ss2.lcl = results$ss2.mean - 1.96*results$ss2.se
  results$ss2.ucl = results$ss2.mean + 1.96*results$ss2.se
  
  results$ss3.mean = mean(squid$ss3)
  results$ss3.se = sd(squid$ss3)
  results$ss3.lcl = results$ss3.mean - 1.96*results$ss3.se
  results$ss3.ucl = results$ss3.mean + 1.96*results$ss3.se
  
  results$rv1.mean = mean(squid$rv1)
  results$rv1.se = sd(squid$rv1)
  results$rv1.lcl = results$rv1.mean - 1.96*results$rv1.se
  results$rv1.ucl = results$rv1.mean + 1.96*results$rv1.se
  
  results$rv2.mean = mean(squid$rv2)
  results$rv2.se = sd(squid$rv2)
  results$rv2.lcl = results$rv2.mean - 1.96*results$rv2.se
  results$rv2.ucl = results$rv2.mean + 1.96*results$rv2.se
  
  results$rv3.mean = mean(squid$rv3)
  results$rv3.se = sd(squid$rv3)
  results$rv3.lcl = results$rv3.mean - 1.96*results$rv3.se
  results$rv3.ucl = results$rv3.mean + 1.96*results$rv3.se
  
  return(results)
}

# finally, a function that iterates through the parameter values
# to estimate mean lambda for each of the 8 population / time combinations

all.models <- function(parameters, it=10000){
  
  for (i in 1:dim(parameters)[1]){
    
    a <- mean.model(parameters[i,], it)
    
    if(i==1){
      output <- a }else{
      
        output <- rbind(output, a)
        
      }
    print(i)
  }
 return(output)
   
}

MPMs.toplot <- all.models(parameters) 

names(MPMs.toplot)[1:3] <- c("stream", "hab", "period")

# save the results
write.csv(MPMs.toplot, file="./output/results/MPMs-toplot.csv", row.names = FALSE)
