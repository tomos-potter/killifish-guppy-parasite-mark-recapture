#=========================================================================================================
# PARASITE PREVALANCE
#=========================================================================================================


data <- H_killi[which(H_killi$sampling>169),] 

data <- as.data.table(data)

data$stream <- as.factor(data$stream)
data$sex <- as.factor(data$sex)
data$Habitat <- as.factor(data$Habitat)


data$count <- 1
data <- data[which(data$Habitat=='C' | data$Habitat=='I'),]

data$habitat2 <- as.character(data$Habitat)

data$size <- ifelse(data$TL>15 & data$TL<=35, "1",
                    ifelse(data$TL>35 & data$TL<=55, "2",
                            ifelse(data$TL>55, "3", NA)))


# infection rate in Caigual
mod.CA <- glmmTMB(parasite ~ -1 +size:habitat2 + (1|sampling), 
                data=data[stream=="CA" & !is.na(size)],
                family="binomial")

predsCA <- emmeans(mod.CA, specs = c("size", "habitat2"), type="response")

# planned contrasts
con <- list(
  size1_cont_int = c(1,0,0,-1,0,0),
  size2_cont_int = c(0,1,0,0,-1,0),
  size3_cont_int = c(0,0,1,0,0,-1)
)

# planned contrast test results:
emmeans(mod.CA, list(~ size:habitat2), contr = con, adjust = "mvt")


killi_ir_CA <- as.data.table(predsCA)
killi_ir_CA$stream <- "CA"


# infection rate in Taylor
mod.TY <- glmmTMB(parasite ~ -1 +size:habitat2 + (1|sampling), 
                  data=data[stream=="TY" & !is.na(size)],
                  family="binomial")


predsTY <- emmeans(mod.TY, specs = c("size", "habitat2"), type="response")
emmeans(mod.TY, list(~ size:habitat2), contr = con, adjust = "mvt")


killi_ir_TY <- as.data.table(predsTY)
killi_ir_TY$stream <- "TY"

# stick them together
killi_ir_focal <- rbind(killi_ir_CA, killi_ir_TY)

killi_ir_focal$streams <- ifelse(killi_ir_focal$stream=="CA", "Caigual", "Taylor")

# save results
write.csv(killi_ir_focal, file="./output/results/killifish-infection-rate-focal.csv", row.names = FALSE)

# infection rate in guppies
gdata <- as.data.table(H_guppy)
gdata <- gdata[(stream=="CA"|stream=="TY") & sampling>169,]

mod.guppy <- glmmTMB(parasite ~ -1 + stream + (1|sampling), 
                    data=gdata,
                    family="binomial")

predsG <- emmeans(mod.guppy, specs = c("stream"), type="response")
guppy_ir <- as.data.table(predsG)
guppy_ir$habitat2 <- "I"
guppy_ir$streams <- ifelse(guppy_ir$stream=="CA", "Caigual", "Taylor")

# save results
write.csv(guppy_ir, file="./output/results/guppy-infection-rate-focal.csv", row.names = FALSE)

# DISSECTION RESULTS
# now lets look at the dissection data

kdata <- dissection_data

kdata$size <- ifelse(kdata$tl>15 & kdata$tl<=35, "1",
                     ifelse(kdata$tl>35 & kdata$tl<=55, "2",
                            ifelse(kdata$tl>55, "3", NA)))

kdata$stream.id <- paste(kdata$drainage, kdata$stream, sep="-")

mod.DIS <- glmmTMB(infected ~ -1+size:community + (1|stream.id), 
                   data=kdata[!is.na(size) & species=="killifish"],
                   family="binomial")

predsDIS <- emmeans(mod.DIS, specs = c("size", "community"), type="response")
emmeans(predsDIS, list(~ size:community), contr = con, adjust = "mvt")


detection_killis <- glmmTMB(parasite_visible_externally ~ tl, 
                            data=kdata[species=="killifish" & !is.na(tl)],
                            family="binomial")

emmeans(detection_killis, specs=c("tl"), type="response")

killi_dis_ir <- as.data.table(predsDIS)
killi_dis_ir$community <- ifelse(killi_dis_ir$community=="KO", "C", "I")

# save results
write.csv(killi_dis_ir, file="./output/results/killifish-infection-rate-dissection.csv", row.names = FALSE)


mod.DIS.g <- glmmTMB(infected ~ 1  + (1|stream.id), 
                     data=kdata[species=="guppy"],
                     family="binomial")

predsDIS.g <- emmeans(mod.DIS.g, specs = "1", type="response")
guppy_dis_ir <- as.data.table(predsDIS.g)

# save results
write.csv(guppy_dis_ir, file="./output/results/guppy-infection-rate-dissection.csv", row.names = FALSE)



