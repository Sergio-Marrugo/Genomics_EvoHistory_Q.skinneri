setwd("/media/biologianeotropical/Sergio_Genomics/Genomics_QS/All_samples/Pop_genomics/PCA")
library(tidyverse)

pca <- read.table("./all_samplesPCA.eigenvec")
eigenval <- scan("./all_samplesPCA.eigenval")
# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# sort out the individual pops
# location
Locality <- rep(NA, length(pca$ind))
Locality[grep("CHA", pca$ind)] <- "Chayotepec"
Locality[grep("COB", pca$ind)] <- "Coban"
Locality[grep("COCH", pca$ind)] <- "Cofradia"
Locality[grep("GDH", pca$ind)] <- "Guevea"
Locality[grep("LAC", pca$ind)] <- "Lashiguxe2"
Locality[grep("LOTO", pca$ind)] <- "Tuxtlas1"
Locality[grep("NIQ", pca$ind)] <- "Niquivil"
Locality[grep("POC", pca$ind)] <- "Pochutla"
Locality[grep("SAC", pca$ind)] <- "San_Antonio"
Locality[grep("SBH", pca$ind)] <- "Honduras"
Locality[grep("SIL", pca$ind)] <- "Lashiguxe1"
Locality[grep("SKC", pca$ind)] <- "Tuxtlas2"
Locality[grep("TPC", pca$ind)] <- "Tres_Picos"
Locality[grep("YAJ", pca$ind)] <- "Yajalon"

locality
# remake data.frame
pca <- as_tibble(data.frame(pca, Locality))
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

#plot pca
b <- ggplot(pca, aes(PC1, PC2, col = Locality)) + geom_point(size = 2)
b <- b + scale_colour_manual(values = c("red", "blue", "orange", "yellow", "purple", "green", "magenta", "cyan", "pink", "black", "cornflowerblue", "deeppink4", "coral3", "darkseagreen" ))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
ggsave("PCA_newlabs.tiff",dpi = 600)


#Only western populations
setwd("D:/Genomics_QS/All_samples/Pop_genomics/PCA")

pca <- read_table2("Western_pops/skinneriPCA_Western.eigenvec", col_names = FALSE)
eigenval <- scan("Western_pops/skinneriPCA_Western.eigenval")
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
Locality <- rep(NA, length(pca$ind))
Locality[grep("CHA", pca$ind)] <- "Chayotepec"
Locality[grep("COCH", pca$ind)] <- "Cofradia"
Locality[grep("GDH", pca$ind)] <- "Guevea"
Locality[grep("LAC", pca$ind)] <- "Lashiguxe2"
Locality[grep("LOTO", pca$ind)] <- "Tuxtlas1"
Locality[grep("SAC", pca$ind)] <- "San_Antonio"
Locality[grep("SIL", pca$ind)] <- "Lashiguxe1"
Locality[grep("SKC", pca$ind)] <- "Tuxtlas2"

# remake data.frame
pca <- as_tibble(data.frame(pca, Locality))

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# plot pca
b <- ggplot(pca, aes(PC1, PC2, col = Locality)) + geom_point(size = 3, shape=19)
b <- b + scale_colour_manual(values = c("red", "orange", "yellow", "green", "magenta", "black", "deeppink4", "coral3"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
ggsave("western_pca.tiff", dpi= 600)