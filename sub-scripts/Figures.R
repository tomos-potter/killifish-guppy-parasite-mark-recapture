#=========================================================================================================
# PLOTTING FIGURES
#=========================================================================================================

#=========================================================================================================
# FIGURE 1 Killifish population densities, survival, and recruitment
#=========================================================================================================

killi <- ggplot(mapping = aes(x = 0:1, y = 1)) +
  theme_void() +
  geom_image(data = data.table(x=0, y=0), 
             aes(image="./data/killifish_silhouette.png", x=x, y=y),
             size=2, inherit.aes = FALSE, colour="darkgrey")

pal <- c("#FFDB6D", "#52854C",
         "#D16103",  "#4E84C4","#C4961A")

toplot <- fread("./output/results/POPAN-parameters.csv")

toplot$streams <- ifelse(toplot$stream=="TY", "Taylor", "Caigual")


toplot$periods <- ifelse(toplot$period==1 & toplot$stream=="CA", "Before",
                         ifelse(toplot$period==2 & toplot$stream=="CA", "Epidemic",
                                ifelse(toplot$period==3 & toplot$stream=="CA", "Endemic",
                                       ifelse(toplot$period==1 & toplot$stream=="TY", "Before 1",
                                              ifelse(toplot$period==2 & toplot$stream=="TY", "Before 2", "Endemic")))))  

toplot$periods<- factor(toplot$periods, levels=c("Before", "Before 1", "Before 2",
                                                 "Epidemic", "Endemic"))                                         

popdens<- ggplot(toplot[parameter=="pop.dens"], aes(x=as.factor(periods), y=log(estimate), colour=habitat, group=habitat))+
  geom_line(linetype="dashed", position = position_dodge(0.25), alpha=0.25) +
  geom_pointrange(aes(ymin=log(lcl), ymax=log(ucl)),
                  position = position_dodge(0.25), size=0.2)+
  theme_tufte()+
  theme(axis.text.x=element_text(size=8), legend.position = "none")+
  geom_text(data = data.frame(streams="Caigual", habitat=NA),
            aes(x=2.5, y=2.75, label="No guppies"), colour="#52854C", size=3)+
  geom_segment(data = data.frame(streams="Caigual", habitat=NA),
               aes(x = 2.1, y = 2.75, xend = 1.5, yend = 2.2),colour="#52854C", size=0.5,arrow = arrow(length = unit(0.3, "cm")))+
  geom_text(data = data.frame(streams="Caigual", habitat=NA),
            aes(x=2, y=-0.75, label="With guppies"), colour="#C4961A", size=3)+
  geom_segment(data = data.frame(streams="Caigual", habitat=NA),
               aes(x = 2, y = -0.5, xend = 1.5, yend = 0.8),colour="#C4961A", size=0.5,arrow = arrow(length = unit(0.3, "cm")))+
  ylim(c(-1,3))+
  labs(x="", y="Log killifish per metre",
       colour="Community", title="A")+
  scale_colour_manual(labels=c("No guppies", "With guppies"),
                      values=pal[c(2,5)])+
  facet_wrap(~streams, scales="free")

surv<- ggplot(toplot[parameter=="Phi"], aes(x=as.factor(periods), y=estimate, colour=habitat, group=habitat))+
  geom_line(linetype="dashed", position = position_dodge(0.25), alpha=0.25) +
  geom_pointrange(aes(ymin=lcl, ymax=ucl),
                  position = position_dodge(0.25), size=0.2)+
  theme_tufte()+
  theme(axis.text.x=element_text(size=8), legend.position="none")+
  ylim(c(0.70,0.95))+
  labs(x="", y="Killifish monthly survival",
       colour="Community", title = "B")+
  scale_colour_manual(labels=c("No guppies", "With guppies"),
                      values=pal[c(2,5)])+
  
  facet_wrap(~streams, scales="free")


rec<- ggplot(toplot[parameter=="babies"], aes(x=as.factor(periods), y=estimate, colour=habitat, group=habitat))+
  geom_line(linetype="dashed", position = position_dodge(0.25), alpha=0.25) +
  geom_pointrange(aes(ymin=lcl, ymax=ucl),
                  position = position_dodge(0.25), size=0.2)+
  theme_tufte()+
  theme(axis.text.x=element_text(size=8), legend.position="none")+
  ylim(c(0.15,0.81))+
  labs(x="Stage of parasite invasion", y="Per-capita monthly recruitment",
       colour="Community", title="C")+
  scale_colour_manual(labels=c("No guppies", "With guppies"),
                      values=pal[c(2,5)])+
  facet_wrap(~streams, scales="free")

killifish_by_period <- killi / popdens / surv / rec / plot_layout(heights = c(0.3,1,1,1))# & theme(legend.position = 'bottom')

ggsave(killifish_by_period,
       filename = "./output/figures/Figure-1.png", 
       width=6, height=8, units = "in")

#=========================================================================================================
# FIGURE 2 - Guppy population densities, survival, and recruitment
#=========================================================================================================

guppy <- ggplot(mapping = aes(x = 0:1, y = 0:1)) +
  theme_void() +
  ylim(c(0,1))+
  xlim(c(0,1))+
  geom_image(data = data.table(x=0.5, y=0.7), 
             aes(image="./data/guppy_silhouette.png", x=x, y=y),
             size=2, inherit.aes = FALSE, colour="darkgrey")

toplot <- fread("./output/results/G-POPAN-parameters.csv")

toplot$streams <- ifelse(toplot$stream=="TY", "Taylor", "Caigual")

toplot$period <- as.factor(ifelse(toplot$period==1 & toplot$stream=="CA", "Before",
                                  ifelse(toplot$period==2 & toplot$stream=="CA", "Epidemic",
                                         ifelse(toplot$period==3 & toplot$stream=="CA","Endemic",
                                                ifelse(toplot$period==1 & toplot$stream=="TY", "Before 1",
                                                       ifelse(toplot$period==2 & toplot$stream=="TY", "Before 2", "Endemic"))))))

toplot$period<- factor(toplot$period, levels=c("Before", "Before 1", "Before 2",
                                               "Epidemic", "Endemic"))

popdens <-ggplot(toplot[parameter=="pop.dens"], aes(x=period, y=log(estimate), group=streams))+
  geom_line(col="grey", linetype="dashed", alpha=0.5)+
  geom_pointrange(aes(ymin=log(lcl), ymax=log(ucl)), size=0.2)+
  labs(title="A", y="Log guppies per metre", 
       x="")+
  ylim(c(1.5,3))+
  theme_tufte()+
  facet_wrap(~streams, scale="free")+
  theme(legend.position="none")

surv <- ggplot(toplot[parameter=="Phi"], aes(x=period, y=estimate, group=stream))+
  geom_line(col="grey", linetype="dashed",  alpha=0.5)+
  geom_pointrange(aes(ymin=lcl, ymax=ucl), size=0.2)+
  theme_tufte()+
  ylim(c(0.7,0.90))+
  labs(title="B", y="Monthly survival", 
       x="")+
  facet_wrap(~streams, scales="free")

recruits <- ggplot(toplot[parameter=="per.capita.recruitment"], aes(x=period, y=estimate, group=streams))+
  geom_line(col="grey", linetype="dashed", alpha=0.5)+
  geom_pointrange(aes(ymin=lcl, ymax=ucl), size=0.2)+
  ylim(c(0.3,0.7))+
  labs(title="C", y="Per capita rectuitment", 
       x="Stage of parasite invasion")+
  theme_tufte()+
  facet_wrap(~streams, scale="free")


guppies_by_period <-  guppy /popdens / surv / recruits /plot_layout(heights = c(0.3,1,1,1))

ggsave(guppies_by_period,
       filename = "./output/figures/Figure-2.png", 
       width=6, height=8, units = "in")

#=========================================================================================================
# FIGURE 3 - size-dependent survival and growth
#=========================================================================================================
toplot<- fread("./output/results/Multistrata-parameters.csv")

toplot$periods <- ifelse(toplot$period==1 & toplot$stream=="CA", "Before",
                         ifelse(toplot$period==2 & toplot$stream=="CA", "Epidemic",
                                ifelse(toplot$period==3 & toplot$stream=="CA", "Endemic",
                                       ifelse(toplot$period==1 & toplot$stream=="TY", "Before 1",
                                              ifelse(toplot$period==2 & toplot$stream=="TY", "Before 2", "Endemic")))))  

toplot$periods<- factor(toplot$periods, levels=c("Before", "Before 1", "Before 2",
                                                 "Epidemic", "Endemic"))

toplot$habitats <- ifelse(toplot$habitat=="C", "No guppies", "With guppies")

toplot$stratum <- ifelse(toplot$stratum=="1to2", "1 to 2", 
                         ifelse(toplot$stratum=="2to3", "2 to 3", toplot$stratum))


pal <- c("#FFDB6D", "#52854C",
         "#D16103",  "#4E84C4","#C4961A")

killi_size <- ggplot(mapping = aes(x = -2:5, y = -1:3)) +
  theme_void() +
  geom_image(data = data.table(x=c(1.5, 2.5, 3.5), y=2), 
             aes(image=c("./data/killifish_silhouette.png",
                         "./data/killifish_silhouette.png",
                         "./data/killifish_silhouette.png"), x=x, y=y),
             size=c(0.5,0.8,1.5), inherit.aes = FALSE, colour="darkgrey")+
  xlim(c(0,5))+
  ylim(c(-1.5,5))+
  annotate("text", x=c(1.5, 2.5, 3.5), y=4.5, label = c("1", "2", "3"))+
  annotate("text", x=c(1.5, 2.5, 3.5), y=0, 
           label = c("15 - 35 mm", "35 - 55 mm", "> 55 mm"), size=3)+
  annotate("text", x=c(1.5, 2.5, 3.5), y=-1, 
           label = c("juvenile", "young adult", "established adult"), size=3)+
  annotate("text", x=0.75, y=c(2), label = c("Size classes"))

label_data <- data.frame(streams="Caigual", 
                         periods = "Before",
                         #   periods=c("Before", "Epidemic", "Endemic"), 
                         habitats=NA)

label_data$periods<- factor(label_data$periods, levels=c("Before", "Epidemic", "Endemic"))

surv_CA <- ggplot(toplot[stream=="CA" & parameter=="S",], 
                  aes(x=as.factor(stratum), y=estimate, colour=habitats, group=habitats))+
  geom_line(linetype="dashed", position = position_dodge(0.25), alpha=0.25)+
  geom_pointrange(aes(ymin=lcl, ymax=ucl),
                  position = position_dodge(0.25), size=0.2)+
  theme_tufte()+
  
  geom_text(data = label_data,
            aes(x=1.5, y=0.5, label="No\nguppies"), colour="#52854C", size=3)+
  geom_segment(data = label_data,
               aes(x = 1.5, y = 0.6, xend = 0.95, yend = 0.88),colour="#52854C", size=0.5,arrow = arrow(length = unit(0.3, "cm")))+
  geom_text(data = label_data,
            aes(x=2.5, y= 0.6, label="With\nguppies"), colour="#C4961A", size=3)+
  geom_segment(data = label_data,
               aes(x = 2.2, y = 0.65, xend = 1.15, yend = 0.92),colour="#C4961A", size=0.5,arrow = arrow(length = unit(0.3, "cm")))+
  
  ylim(c(0.3, 1))+
  labs(x="", y="Survival",
       colour="Community", title="A", subtitle = "Caigual")+
  scale_colour_manual(values=pal[c(2,5)])+
  facet_wrap(~periods, scales="free")


surv_TY <- ggplot(toplot[stream=="TY" & parameter=="S",], 
                  aes(x=as.factor(stratum), y=estimate, colour=habitats, group=habitats))+
  geom_line(linetype="dashed", position = position_dodge(0.25), alpha=0.25)+
  geom_pointrange(aes(ymin=lcl, ymax=ucl),
                  position = position_dodge(0.25), size=0.2)+
  theme_tufte()+
  ylim(c(0.3, 1))+
  labs(x="Size class", y="Survival",
       colour="Community", title="B", subtitle = "Taylor")+
  scale_colour_manual(labels=c("No guppies", "With guppies"),
                      values=pal[c(2,5)])+
  facet_wrap(~periods, scales="free")

growth_CA <- ggplot(toplot[stream=="CA" & parameter=="Psi" & estimate!=0,], 
                    aes(x=as.factor(stratum), y=estimate, colour=habitats, group=habitats))+
  geom_line(linetype="dashed", position = position_dodge(0.25), alpha=0.25)+
  geom_pointrange(aes(ymin=lcl, ymax=ucl),
                  position = position_dodge(0.25), size=0.2)+
  ylim(c(0,0.8))+
  theme_tufte()+
  labs(x="", y="Growth probability",
       colour="Community", title="C", subtitle = "Caigual")+
  scale_colour_manual(labels=c("No guppies", "With guppies"),
                      values=pal[c(2,5)])+
  facet_wrap(~periods, scales="free")

growth_TY <- ggplot(toplot[stream=="TY" & parameter=="Psi" & estimate!=0,], 
                    aes(x=as.factor(stratum), y=estimate, colour=habitats, group=habitats))+
  geom_line(linetype="dashed", position = position_dodge(0.25), alpha=0.25)+
  geom_pointrange(aes(ymin=lcl, ymax=ucl),
                  position = position_dodge(0.25), size=0.2)+
  ylim(c(0,0.8))+
  theme_tufte()+
  labs(x="Size class transition", y="Growth probability",
       colour="Community", title="D", subtitle = "Taylor")+
  scale_colour_manual(labels=c("No guppies", "With guppies"),
                      values=pal[c(2,5)])+
  facet_wrap(~periods, scales="free")

size_results <-(       killi_size /
                         plot_spacer() /
                         surv_CA / surv_TY / 
                         growth_CA / growth_TY / 
                         plot_layout(heights = c(0.5, 0.1, 1,  1,1,1)) 
                       & theme(legend.position = 'none')) 

ggsave(size_results,
       filename = "./output/figures/Figure-3.png", 
       width=8, height=12, units = "in")


#=========================================================================================================
# FIGURE 4 - infection rates from focal streams and dissection study
#=========================================================================================================

killi_ir_focal <- fread("./output/results/killifish-infection-rate-focal.csv")
killi_dis_ir <- fread("./output/results/killifish-infection-rate-dissection.csv")

# pretty colours
pal <- c("#FFDB6D", "#52854C",
         "#D16103",  "#4E84C4","#C4961A")


parasite.plot <- ggplot(killi_ir_focal, aes(x=as.factor(size), y=prob, group=habitat2, colour=habitat2))+
  geom_text(data = data.frame(streams="Caigual", habitat2=NA),
            aes(x=1, y=0.1, label="No\nguppies"), colour="#52854C", size=3)+
  geom_segment(data = data.frame(streams="Caigual", habitat2=NA),
               aes(x = 1.15, y = 0.1, xend = 1.5, yend = 0.15),colour="#52854C", size=0.5,arrow = arrow(length = unit(0.3, "cm")))+
  geom_text(data = data.frame(streams="Caigual", habitat2=NA),
            aes(x=2.5, y=0.5, label="With\nguppies"), colour="#C4961A", size=3)+
  geom_segment(data = data.frame(streams="Caigual", habitat2=NA),
               aes(x = 2.5, y = 0.4, xend = 2.5, yend = 0.2),colour="#C4961A", size=0.5,arrow = arrow(length = unit(0.3, "cm")))+
  geom_line(linetype="dashed", position = position_dodge(0.25), alpha=0.25)+
  geom_pointrange(aes(ymin=asymp.LCL, ymax=asymp.UCL),
                  position = position_dodge(0.25), size=0.2)+
  theme_tufte()+
  labs(x="Size class", y="Percent visibly parasitized",
       colour="Community", title="A")+
  scale_x_discrete(labels=c("1", "2","3" ))+
  scale_y_continuous(labels = scales::percent, limits = c(0,1))+
  scale_colour_manual(labels=c("No guppies", "With guppies"),
                      values=pal[c(2,5)])+
  facet_wrap(~streams)



dis.plot <- ggplot(killi_dis_ir, aes(x=as.factor(size), y=prob, group=community, colour=community))+
  geom_line(linetype="dashed", position = position_dodge(0.25), alpha=0.25)+
  theme_tufte()+
  geom_pointrange(aes(ymin=asymp.LCL, ymax=asymp.UCL),
                  position = position_dodge(0.25), size=0.2)+
  labs(x="Size class", y="Percent parasitized on dissection",
       colour="Community", title="B")+
  scale_x_discrete(labels=c("1", "2","3" ))+
  scale_y_continuous(labels = scales::percent, limits = c(0,1))+
  scale_colour_manual(labels=c("No guppies", "With guppies"),
                      values=pal[c(2,5)])


comb.parasite.plot <- parasite.plot + dis.plot +
  plot_layout(guides="collect", widths = c(2,1)) & theme(legend.position = 'none')

ggsave(comb.parasite.plot,
       filename = "./output/figures/Figure-4.png", 
       width=8, height=3.5, units = "in")

#=========================================================================================================
# FIGURE 5 Killifish population growth rates (MPM predicted lambdas)
#=========================================================================================================

MPMs.toplot <- fread("./output/results/MPMs-toplot.csv")

MPMs.toplot$periods <- ifelse(MPMs.toplot$period==1 & MPMs.toplot$stream=="CA", "Before",
                              ifelse(MPMs.toplot$period==2 & MPMs.toplot$stream=="CA", "Epidemic",
                                     ifelse(MPMs.toplot$period==3 & MPMs.toplot$stream=="CA", "Endemic",
                                            ifelse(MPMs.toplot$period==1 & MPMs.toplot$stream=="TY", "Before 1",
                                                   ifelse(MPMs.toplot$period==2 & MPMs.toplot$stream=="TY", "Before 2", "Endemic")))))  

MPMs.toplot$periods<- factor(MPMs.toplot$periods, levels=c("Before", "Before 1", "Before 2",
                                                           "Epidemic", "Endemic"))     


MPMs.toplot$habitats <- ifelse(MPMs.toplot$hab=="C", "No guppies", "With guppies")
MPMs.toplot$streams <- ifelse(MPMs.toplot$stream=="TY", "Taylor", "Caigual")

pal <- c("#FFDB6D", "#52854C",
         "#D16103",  "#4E84C4","#C4961A")

MPMfig <- ggplot(MPMs.toplot, aes(x=periods, y=lambda.mean, colour=habitats))+
  geom_hline(yintercept = 1, linetype = "dotted")+
  geom_pointrange(aes(ymin=lambda.lcl, ymax=lambda.ucl),
                  position = position_dodge(0.5))+
  geom_text(data = data.frame(streams="Caigual", habitats=NA),
            aes(x=1, y=0.95, label="No\nguppies"), colour="#52854C", size=3)+
  geom_segment(data = data.frame(streams="Caigual", habitats=NA),
               aes(x = 1, y = 0.975, xend = 0.9, yend = 1.05),colour="#52854C", size=0.5,arrow = arrow(length = unit(0.3, "cm")))+
  geom_text(data = data.frame(streams="Caigual", habitats=NA),
            aes(x=1.75, y=0.9, label="With\nguppies"), colour="#C4961A", size=3)+
  geom_segment(data = data.frame(streams="Caigual", habitats=NA),
               aes(x = 1.75, y = 0.925, xend = 1.25, yend = 1.05),colour="#C4961A", size=0.5,arrow = arrow(length = unit(0.3, "cm")))+
  theme_tufte()+          
  theme(legend.position = "none")+
  scale_colour_manual(values=pal[c(2,5)])+
  ylim(c(0.86, 1.14))+
  labs(x="Stage of parasite invasion", y="Population growth rate λ",
       colour="Community")+
  facet_wrap(~streams, scales="free")

ggsave(MPMfig,        
       filename = "./output/figures/Figure-5.png", 
       width=8, height=3, units = "in")

# ANNOUNCE THAT WE HAVE REACHED THE END OF THE CODE!
system("say stick a fork in me I am done")
