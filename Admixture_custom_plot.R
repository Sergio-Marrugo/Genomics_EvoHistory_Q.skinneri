library(ggplot2)
library(reshape2)
library(stringr)
library(dplyr)
setwd("C:/Users/EQUIPO/Desktop/Cuarto_semestre/Documento_Tesis/admixture_map")

#load the ind_pop map
inds_pops <- read.table("inds_pops.txt")
colnames(inds_pops) <- c("individuo", "poblacion")

#Load the K=5 admiture results
admixture<- read.table("K5_Q_to_plot")
colnames(admixture)<- c("individuo", "c1", "c2", "c3", "c4", "c5")
admixture <- mutate(admixture, c6 = 0)
admixture <- admixture %>% relocate(c6, .before = c5)
admix_long <- inner_join(admixture,inds_pops, by = "individuo")
#change group CHM by SNO
admix_long$poblacion <- str_replace(admix_long$poblacion, "CHM", "SNO")
write.csv(admix_long, file = "admix_plot_table_k5.csv")
#file to plot
arranged_k5<-read.csv("admix_plot_table_k5.csv")




#test 
d.test <- data.frame(A1=c(0.000, 0.000, 0.020, 0.000, 0.000, 0.000),
                     A2=c(0.000, 0.000, 0.235, 0.195, 0.166, 0.205),
                     A3=c(0.065, 0.027, 0.000, 0.027, 0.000, 0.036),
                     A4=c(0.000, 0.000, 0.007, 0.011, 0.000, 0.000),
                     A5=c(0.000, 0.000, 0.000, 0.002, 0.028, 0.000),
                     A6=c(0.000, 0.041, 0.021, 0.068, 0.106, 0.105),
                     A7=c(0.093, 0.085, 0.001, 0.056, 0.110, 0.000),
                     A8=c(0.000, 0.000, 0.000, 0.000, 0.000, 0.029),
                     A9=c(0.000, 0.000, 0.058, 0.027, 0.096, 0.156),
                     A10=c(0.000, 0.023, 0.129, 0.012, 0.074, 0.117),
                     A11=c(0.000, 0.041, 0.000, 0.000, 0.000, 0.000),
                     A12=c(0.024, 0.000, 0.000, 0.000, 0.000, 0.000),
                     A13=c(0.817, 0.783, 0.527, 0.446, 0.258, 0.321),
                     A14=c(0.000, 0.000, 0.000, 0.006, 0.000, 0.000),
                     A15=c(0.000, 0.000, 0.000, 0.054, 0.143, 0.027),
                     A16=c(0.000, 0.000, 0.000, 0.000, 0.000, 0.003),
                     A17=c(0.000, 0.000, 0.000, 0.097, 0.019, 0.000))

# Add id and population label columns. Needed for melting and plotting.
d.test$population = c("p1", "p1", "p2", "p2", "p2", "p2")
d.test$subject_id = paste("id", 1:6, sep="")

# Melt (reshape data from wide format to long format).
mdat = melt(d.test, id.vars=c("subject_id", "population"), 
            variable.name="Ancestry", value.name="Fraction")
#let's try with our data

extra_long_admix <- melt(arranged_k5, id.vars=c("individuo", "poblacion"), 
                         variable.name="Ancestry", value.name="Fraction")

#let's plot that shiiiii
ggplot(extra_long_admix, aes(x=individuo, y=Fraction, fill=Ancestry)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(. ~ poblacion, drop=TRUE, space="free", scales="free")+
  theme(strip.background=element_blank())
  #theme(strip.text=element_text(size=12))

ggsave("K5_admix_plot.jpg")

#let's plot k=6
k6_table <- read.table("skinneriAdmix.6.Q")
colnames(k6_table)<- c("c1", "c2", "c3", "c4", "c5", "c6")
K6_admix <- cbind(inds_pops,k6_table)
write.csv(K6_admix, "admix_plot_table_k6.csv")
arranged_k6 <- read.csv("admix_plot_table_k6.csv")

k6_long_admix <- melt(arranged_k6, id.vars=c("individuo", "poblacion"), 
                         variable.name="Ancestry", value.name="Fraction")

ggplot(k6_long_admix, aes(x=individuo, y=Fraction, fill=Ancestry)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(. ~ poblacion, drop=TRUE, space="free", scales="free")+
  theme(strip.background=element_blank())

