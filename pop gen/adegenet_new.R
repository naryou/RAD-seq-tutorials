getwd()
setwd('/Users/narcis/Narcis/Project/Courses_Conferences/Bioinformatics_NY/PartII/R_course/')

#install.packages('adegenet', dep=TRUE) #dependencies TRUE
#install.packages('vcfR')
#install.packages('devtools')
#devtools::install_github("thierrygosselin/radiator") # I installed radiator to use its genomic_convertor function

library(adegenet)
library(ape)
library(pegas) #this might not be needed for PCA
library(seqinr) #this might not be needed for PCA
library(ggplot2) #this might not be needed for PCA
library(vcfR)
#library(help = vcfR)
#library(radiator)

#packageDescription("adegenet", fields = "Version")
??adegenet
?df2genind
#genomic_converter() # this is to make genlight object but vcfR does it.
#df2genind(X, sep = NULL, ncode = NULL, ind.names = NULL,
#          loc.names = NULL, pop = NULL, NA.char = "", ploidy = 2,
#          type = c("codom", "PA"), strata = NULL, hierarchy = NULL)


# df <- read.vcfR('M:/Project/R/Rworkspace_M/adegenet/warn/m3M0n1N0_80/p1r5maf05.vcf', limit = 1e+07, nrows = -1, skip = 0, cols = NULL, convertNA = TRUE, verbose = TRUE)

#df2 <- read.table('./variants_miss100.recode.vcf', header = T, sep = '\t')
df <- read.vcfR('./populations.snps.vcf', limit =20000000, nrows = -1, skip = 0, cols = NULL, convertNA = TRUE, verbose = TRUE)
df_genind <- vcfR2genind(df, sep = "[/]") 


#NB: another way to read data in is import2genind() from adegenet. it tries to detect the file format and finds the appropriate import function
# for import2genind(), structure should be with .str or .stru, and genepop with .gen extention
# you can also use read.structure or read.genepop

#PCA
X <- scaleGen(df_genind, NA.method="mean") #Replaces missing data with mean of allele frequencies
class(X)
dim(X)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3) #dudi.pco performs PCoA, it needs a distance matrix
#
s.label(pca1$li)
title("PCA of GSP_47")
add.scatter.eig(pca1$eig[1:20], 3,1,2)
#
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50)) #Good!
tiff("./p1r5maf05_2.tiff", width = 7.5, height = 3.5, pointsize = 1/1000, units = 'in', res = 1000)
par(mgp=c(3,0.8,0),mar=c(4,4,1,1)+0.1)
plot(pca1$li[,1], pca1$li[,2] ,  cex=2,lwd=3, pch = c(2,1,3,0), #col=color_ny, # #$li is the principal components of the analysis; 
     cex.lab=1, cex.axis=1, 
     ylab=paste("PC2", paste("(", round(pca1$eig[2]/sum(pca1$eig)*100), " %)", sep="" )),
     xlab=paste("PC1", paste("(", round(pca1$eig[1]/sum(pca1$eig)*100), " %)", sep="" )))
#text(pca1$li[,1], pca1$li[,2], labels=warn_info$Country, cex= 0.4) # YAAAYYYY
legend("bottomright",legend=paste(c('E_Canada', 'Europe', 'Russia', 'W_Canada')),
       pch = c(2,1,3,0), cex=0.6)

dev.off()
#
#
#
df_genlight <- vcfR2genlight(df, n.cores = 2) 
??vcfR2genlight
??glPca
??glDotProd
pca1 <- glPca(df_genlight)
pca1 <- glPca(df_genind)
df <- read.vcfR('./variants_miss75_SNP_only.recode.vcf', limit =20000000, nrows = -1, skip = 0, cols = NULL, convertNA = TRUE, verbose = TRUE)
# df_genind <- vcfR2genind(df, sep = "[/]") 

# everything related to genlight class is following http://adegenet.r-forge.r-project.org/files/tutorial-genomics.pdf
df_genlight <- vcfR2genlight(df, parallel = TRUE) 
#
pca1 <- glPca(df_genlight)
pca1
scatter(pca1, posi="bottomright") 
title("PCA of the rubellum\n axes 1-2")
# to add color
myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4) 
abline(h=0,v=0, col="grey") 
add.scatter.eig(pca1$eig[1:40],2,1,2, posi="topright", inset=.05, ratio=.3)

# to confirm the pca with NJ tree:
library(ape) 
tre <- nj(dist(as.matrix(df_genind))) 
tre

plot(tre, typ="fan", cex=0.7) 
title("NJ tree of the rubellum data")

# tree with color
plot(tre, typ="fan", show.tip=FALSE) 
tiplabels(pch=20, col=myCol, cex=4) 
title("NJ tree of the rubellum data")

# DAPC
dapc1 <- dapc(df_genind, n.pca=10, n.da=1) # let's see if it works
scatter(dapc1,scree.da=FALSE, bg="white", posi.pca="topright", legend=TRUE, 
        txt.leg=paste("group", 1:2), col=c("red","blue"))
compoplot(dapc1, col=c("red","blue"),lab="", txt.leg=paste("group", 1:2), ncol=2)
loadingplot(dapc1$var.contr, thres=1e-3)
loadingplot(tail(dapc1$var.contr[,1],100), thres=1e-3) #this is not possible with 1M SNPs!!!

#
#
rub_info <- read.table('./rubxxxxxx.txt', header = T, sep = '\t')
colnames(rub_info)
#
col <- rub_info$Cluster
rbPal <- colorRampPalette(c('red3','green', 'blue3')) # 'blue3','green','red3'
color_ny <- rbPal(20)[as.numeric(cut(col,breaks = 20))]
shapes <- rub_info$Region_L
levels(shapes)
#
#autoplot(prcomp(X), data= gcf_clusters, colour = 'Cluster') ##load ggfortify for this

#PCA
X <- scaleGen(df_genlight, NA.method="mean") #Replaces missing data with mean of allele frequencies
class(X)
dim(X)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3) #dudi.pco performs PCoA, it needs a distance matrix
#
s.label(pca1$li)
title("PCA of GSP_47")
add.scatter.eig(pca1$eig[1:20], 3,1,2)
#
#barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50)) #Good!
tiff("M:/Project/R/Rworkspace_M/adegenet/GSP/fallax_m3_pop_47/p1r5maf05_2.tiff", width = 7.5, height = 3.5, pointsize = 1/1000, units = 'in', res = 1000)
par(mgp=c(3,0.8,0),mar=c(4,4,1,1)+0.1)
plot(pca1$li[,1], pca1$li[,2] ,  cex=2,lwd=3, col=color_ny , pch = c(2,1,3,0)[as.numeric(shapes)], #col=color_ny, # #$li is the principal components of the analysis; 
     cex.lab=1, cex.axis=1, 
     ylab=paste("PC2", paste("(", round(pca1$eig[2]/sum(pca1$eig)*100), " %)", sep="" )),
     xlab=paste("PC1", paste("(", round(pca1$eig[1]/sum(pca1$eig)*100), " %)", sep="" )))
#text(pca1$li[,1], pca1$li[,2], labels=warn_info$Country, cex= 0.4) # YAAAYYYY
legend("bottomright",legend=paste(c('E_Canada', 'Europe', 'Russia', 'W_Canada')),
       pch = c(2,1,3,0), cex=0.6)

dev.off()
#

#
mag_info <- read.table('./adegenet/GSP/p1r5maf05_randSNP/all_info.txt', header = T, sep = '\t')
colnames(mag_info)
df <- read.vcfR('M:/Project/R/Rworkspace_M/adegenet/GSP/p1r5maf05_randSNP/batch_1.vcf', limit = 1e+07, nrows = -1, skip = 0, cols = NULL, convertNA = TRUE, verbose = TRUE)
# limit: amount of memory (in bytes) not to exceed when reading in a ???le
# ###### for GSP_all ------------------------------------------------------

df_genind <- vcfR2genind(df, sep = "[/]") 
#
col <- mag_info$Cluster
rbPal <- colorRampPalette(c('red3','green', 'blue3')) # 'blue3','green','red3'
color_ny <- rbPal(20)[as.numeric(cut(col,breaks = 20))]
shapes <- mag_info$Region_L
levels(shapes)
#
#autoplot(prcomp(X), data= gcf_clusters, colour = 'Cluster') ##load ggfortify for this

#PCA
X <- scaleGen(df_genind, NA.method="mean") #Replaces missing data with mean of allele frequencies
class(X)
dim(X)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3) #dudi.pco performs PCoA, it needs a distance matrix
#
s.label(pca1$li)
title("PCA of GSP_47")
add.scatter.eig(pca1$eig[1:20], 3,1,2)
#
#barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50)) #Good!
tiff("M:/Project/R/Rworkspace_M/adegenet/GSP/fallax_m3_pop_47/p1r5maf05_2.tiff", width = 7.5, height = 3.5, pointsize = 1/1000, units = 'in', res = 1000)
par(mgp=c(3,0.8,0),mar=c(4,4,1,1)+0.1)
plot(pca1$li[,1], pca1$li[,2] ,  cex=2,lwd=3, col=color_ny , pch = c(2,1,3,0)[as.numeric(shapes)], #col=color_ny, # #$li is the principal components of the analysis; 
     cex.lab=1, cex.axis=1, 
     ylab=paste("PC2", paste("(", round(pca1$eig[2]/sum(pca1$eig)*100), " %)", sep="" )),
     xlab=paste("PC1", paste("(", round(pca1$eig[1]/sum(pca1$eig)*100), " %)", sep="" )))
#text(pca1$li[,1], pca1$li[,2], labels=warn_info$Country, cex= 0.4) # YAAAYYYY
legend("bottomright",legend=paste(c('E_Canada', 'Europe', 'Russia', 'W_Canada')),
       pch = c(2,1,3,0), cex=0.6)

dev.off()
#
#
####### for GSP_47
mag_info <- read.table('./adegenet/GSP/fallax_m3_pop_47/info_47.txt', header = T, sep = '\t')
colnames(mag_info)
df <- read.vcfR('M:/Project/R/Rworkspace_M/adegenet/GSP/fallax_m3_pop_47/p1r5maf05_47.vcf', limit = 1e+07, nrows = -1, skip = 0, cols = NULL, convertNA = TRUE, verbose = TRUE)

df_genind <- vcfR2genind(df, sep = "[/]") 
#
col <- mag_info$Cluster
rbPal <- colorRampPalette(c('red3','green', 'blue3')) # 'blue3','green','red3'
color_ny <- rbPal(20)[as.numeric(cut(col,breaks = 20))]
shapes <- mag_info$Region_L
levels(shapes)
#
#autoplot(prcomp(X), data= gcf_clusters, colour = 'Cluster') ##load ggfortify for this

#PCA
X <- scaleGen(df_genind, NA.method="mean") #Replaces missing data with mean of allele frequencies
class(X)
dim(X)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3) #dudi.pco performs PCoA, it needs a distance matrix
#
s.label(pca1$li)
title("PCA of GSP_47")
add.scatter.eig(pca1$eig[1:20], 3,1,2)
#
#barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50)) #Good!
tiff("M:/Project/R/Rworkspace_M/adegenet/GSP/fallax_m3_pop_47/p1r5maf05_2.tiff", width = 7.5, height = 3.5, pointsize = 1/1000, units = 'in', res = 1000)
par(mgp=c(3,0.8,0),mar=c(4,4,1,1)+0.1)
plot(pca1$li[,1], pca1$li[,2] ,  cex=2,lwd=3, col=color_ny , pch = c(2,1,3,0)[as.numeric(shapes)], #col=color_ny, # #$li is the principal components of the analysis; 
     cex.lab=1, cex.axis=1, 
     ylab=paste("PC2", paste("(", round(pca1$eig[2]/sum(pca1$eig)*100), " %)", sep="" )),
     xlab=paste("PC1", paste("(", round(pca1$eig[1]/sum(pca1$eig)*100), " %)", sep="" )))
#text(pca1$li[,1], pca1$li[,2], labels=warn_info$Country, cex= 0.4) # YAAAYYYY
legend("bottomright",legend=paste(c('E_Canada', 'Europe', 'Russia', 'W_Canada')),
 pch = c(2,1,3,0), cex=0.6)

dev.off()
#
##DAPC
grp <- find.clusters(ref, max.n.clust=10, n.pca=200)
3
str(grp)
is(grp)
length(grp[[3]])
dapc1 <- dapc(ref, grp$grp, n.pca = 3, n.da=2) #group by cluster
scatter(dapc1)
#
#
#

# #For global/association data --------------------------------------------
mag_info <- read.table('./adegenet/GSP/fallax_m3_exma/info_exma.txt', header = T, sep = '\t')
colnames(mag_info)
df <- read.vcfR('M:/Project/R/Rworkspace_M/adegenet/GSP/fallax_m3_exma/p1r5maf05.vcf', limit = 1e+07, nrows = -1, skip = 0, cols = NULL, convertNA = TRUE, verbose = TRUE)

df_genind <- vcfR2genind(df, sep = "[/]") 
#
col <- mag_info$Genetic_cluster
rbPal <- colorRampPalette(c('blue3','green','red3')) # 'blue3','green','red3'
color_ny <- rbPal(20)[as.numeric(cut(col,breaks = 20))]
shapes <- mag_info$Region_L
levels(shapes)
#
#autoplot(prcomp(X), data= gcf_clusters, colour = 'Cluster') ##load ggfortify for this

#PCA
X <- scaleGen(df_genind, NA.method="mean") #Replaces missing data with mean of allele frequencies
class(X)
dim(X)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3) #dudi.pco performs PCoA, it needs a distance matrix
#
s.label(pca1$li)
title("PCA of GSP")
add.scatter.eig(pca1$eig[1:20], 3,1,2)
#
#barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50)) #Good!
tiff("M:/Project/R/Rworkspace_M/adegenet/GSP/fallax_m3_exma/p1r5maf05.tiff", width = 6, height = 3, pointsize = 1/1000, units = 'in', res = 1000)
par(mgp=c(3,0.8,0),mar=c(4,4,1,1)+0.1)
plot(pca1$li[,1], pca1$li[,2] ,  cex=2,lwd=3, col=color_ny , pch = c(2,1,3,0)[as.numeric(shapes)], #col=color_ny, # #$li is the principal components of the analysis; 
     cex.lab=1, cex.axis=1, ylab=paste("PC2", paste("(", round(pca1$eig[2]/sum(pca1$eig)*100), " %)", sep="" )),
     xlab=paste("PC1", paste("(", round(pca1$eig[1]/sum(pca1$eig)*100), " %)", sep="" )))
#text(pca1$li[,1], pca1$li[,2], labels=warn_info$Country, cex= 0.4) # YAAAYYYY
#legend("bottomleft",legend=paste(c('Arctic', 'Boreal', 'Temperate')),
     #  pch = c(2,1,3), cex=0.6)

dev.off()
#
##DAPC
grp <- find.clusters(ref, max.n.clust=10, n.pca=200)
3
str(grp)
is(grp)
length(grp[[3]])
dapc1 <- dapc(ref, grp$grp, n.pca = 3, n.da=2) #group by cluster
scatter(dapc1)
#
#
#

# ANNOVAR_sig_157 SNPs ----------------------------------------------------
mag_info <- read.table('./adegenet/GSP/fallax_m3_exma/info_exma.txt', header = T, sep = '\t')
colnames(mag_info)
#
col <- mag_info$Genetic_cluster
rbPal <- colorRampPalette(c('blue3','green','red3')) # 'blue3','green','red3'
color_ny <- rbPal(20)[as.numeric(cut(col,breaks = 20))]
shapes <- mag_info$Region_L
levels(shapes)

df <- read.vcfR('M:/Project/R/Rworkspace_M/adegenet/GSP/ANNOVAR_sig_157/can_SNPs_all_ld.3_cor_rm_157.vcf', limit = 1e+07, nrows = -1, skip = 0, cols = NULL, convertNA = TRUE, verbose = TRUE)

df_genind <- vcfR2genind(df, sep = "[/]") 

#PCA
X <- scaleGen(df_genind, NA.method="mean") #Replaces missing data with mean of allele frequencies
class(X)
dim(X)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3) #dudi.pco performs PCoA, it needs a distance matrix
#
s.label(pca1$li)
title("PCA of GSP_157 SNPs")
add.scatter.eig(pca1$eig[1:20], 3,1,2)
#
#barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50)) #Good!
??tiff
tiff("M:/Project/R/Rworkspace_M/adegenet/GSP/fallax_m3_exma/p1r5maf05.tiff", width = 6, height = 3, pointsize = 1/1000, units = 'in', res = 1000)
par(mgp=c(3,0.8,0),mar=c(4,4,1,1)+0.1)
plot(pca1$li[,1], pca1$li[,2] ,  cex=2,lwd=3, col=color_ny , pch = c(2,1,3,0)[as.numeric(shapes)], #col=color_ny, # #$li is the principal components of the analysis; 
     cex.lab=1, cex.axis=1, ylab=paste("PC2", paste("(", round(pca1$eig[2]/sum(pca1$eig)*100), " %)", sep="" )),
     xlab=paste("PC1", paste("(", round(pca1$eig[1]/sum(pca1$eig)*100), " %)", sep="" )))
legend("bottomleft",legend=paste(c("E_Canada","Europe","Russia","W_Canada")),
  pch = c(2,1,3,0), cex=0.6)

dev.off()
#
##DAPC
grp <- find.clusters(ref, max.n.clust=10, n.pca=200)
3
str(grp)
is(grp)
length(grp[[3]])
dapc1 <- dapc(ref, grp$grp, n.pca = 3, n.da=2) #group by cluster
scatter(dapc1)
#
#

###########
#PCA
X <- scaleGen(ref, NA.method="mean") #Replaces missing data with mean of allele frequencies
class(X)
dim(X)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3) #dudi.pco performs PCoA, it needs a distance matrix
#
s.label(pca1$li)
title("PCA of GSP magellanicum")
add.scatter.eig(pca1$eig[1:20], 3,1,2)
#
#s.class(pca1$li, pop(microbov)) #If I had pop info, I could use this
#title("PCA of microbov datasetnnaxes 1-2")
#add.scatter.eig(pca1$eig[1:20], 3,1,2)
#
#barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50)) #Good!
plot(pca1$li[,1:2], col=color_ny, pch = c(0, 2, 3, 1, 8)[as.numeric(shapes)], cex=3,lwd=5, #$li is the principal components of the analysis;
     ylab=paste("Axis 2", paste("(", round(pca1$eig[2]/sum(pca1$eig)*100), " %)", sep="" )), #I don't understand the axis
     xlab=paste("Axis 1", paste("(", round(pca1$eig[1]/sum(pca1$eig)*100), " %)", sep="" )))

legend("topleft",legend=paste(c("E_Canada","Europe","Russia","Scandinavia", "W_Canada")),
       pch = c(0, 2, 3, 1, 8))

levels(mag_coord$Region_L)
col <- mag_coord$Genetic_cluster #It should be in the same order as the genepop file
rbPal <- colorRampPalette(c('blue3','red3'))
color_ny <- rbPal(20)[as.numeric(cut(col,breaks = 20))]
shapes <- mag_coord$Region_L
#

#ggplot(pca1, aes(x=eig[1], y=eig[2])) + geom_point(size=3) #Argggg, it doesn't work


# #global_all ------GCF & GSP combined-------------------------------------------------------

####### for global_all
mag_info <- read.table('./adegenet/GSP/p1r5maf05_randSNP/xxxxx_info.txt', header = T, sep = '\t')
colnames(mag_info)
df <- read.vcfR('M:/Project/R/Rworkspace_M/adegenet/global_all/ref_m3_all_p1r5maf05_rand/batch_1.vcf', limit = 1e+07, nrows = -1, skip = 0, cols = NULL, convertNA = TRUE, verbose = TRUE)

df_genind <- vcfR2genind(df, sep = "[/]") 
#
col <- mag_info$Cluster
rbPal <- colorRampPalette(c('red3','green', 'blue3')) # 'blue3','green','red3'
color_ny <- rbPal(20)[as.numeric(cut(col,breaks = 20))]
shapes <- mag_info$Region_L
levels(shapes)
#
#autoplot(prcomp(X), data= gcf_clusters, colour = 'Cluster') ##load ggfortify for this

#PCA
X <- scaleGen(df_genind, NA.method="mean") #Replaces missing data with mean of allele frequencies
class(X)
dim(X)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3) #dudi.pco performs PCoA, it needs a distance matrix
#
s.label(pca1$li)
title("PCA of global_all")
add.scatter.eig(pca1$eig[1:20], 3,1,2)
#
#barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50)) #Good!
tiff("M:/Project/R/Rworkspace_M/adegenet/GSP/fallax_m3_pop_47/p1r5maf05_2.tiff", width = 7.5, height = 3.5, pointsize = 1/1000, units = 'in', res = 1000)
par(mgp=c(3,0.8,0),mar=c(4,4,1,1)+0.1)
plot(pca1$li[,1], pca1$li[,2] ,  cex=2,lwd=3, col=color_ny , pch = c(2,1,3,0)[as.numeric(shapes)], #col=color_ny, # #$li is the principal components of the analysis; 
     cex.lab=1, cex.axis=1, 
     ylab=paste("PC2", paste("(", round(pca1$eig[2]/sum(pca1$eig)*100), " %)", sep="" )),
     xlab=paste("PC1", paste("(", round(pca1$eig[1]/sum(pca1$eig)*100), " %)", sep="" )))
#text(pca1$li[,1], pca1$li[,2], labels=warn_info$Country, cex= 0.4) # YAAAYYYY
legend("bottomright",legend=paste(c('E_Canada', 'Europe', 'Russia', 'W_Canada')),
       pch = c(2,1,3,0), cex=0.6)

dev.off()



# #For local data/paper -----------------------------------------------------
df_local <- read.table('./adegenet/local/paper_local_p1r5maf05.genepop', header = F, sep = '\t') #modified genepop file > no pop, no header
Clusters <- read.table('./adegenet/local/clusters.txt', header = T, sep = '\t')
df_local$V1 <- NULL #To get rid of the ind names. it seems like the software reads them as the alleles
df_local <- apply(df_local, 2, function(x) as.character(x))
df_local[5:5]
df_local[df_local=="0"] <- NA

df_local_gen <- df2genind(df_local, ploidy = 1, ncode = 3)
head(df_local_gen)
head(df1_gen) #Good, but why it says 1 ind?
#tab(df1_gen)
#sum(is.na(df1_gen$tab))
##DAPC
grp <- find.clusters(df_local_gen, max.n.clust=10, n.pca=200) #NY: this was grp_NY in the code, I changed it to grp, bcs the next is grp...
2
str(grp)
is(grp)
length(grp[[3]])
dapc1 <- dapc(df_local_gen, grp$grp, n.pca = 3, n.da=2) #group by cluster
scatter(dapc1)
#PCA
X <- scaleGen(df_local_gen, NA.method="mean") #Replaces missing data with mean of allele frequencies
class(X)
dim(X)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3) #dudi.pco performs PCoA, it needs a distance matrix
#
s.label(pca1$li)
title("PCA of local magellanicum")
add.scatter.eig(pca1$eig[1:20], 3,1,2)
#
#s.class(pca1$li, pop(microbov)) #If I had pop info, I could use this
#title("PCA of microbov datasetnnaxes 1-2")
#add.scatter.eig(pca1$eig[1:20], 3,1,2)
#
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50)) #Good!
tiff("./adegenet/local/local_adegenet_6-3.5.tiff", width = 6, height = 3.5, pointsize = 1/1000, units = 'in', res = 1000)
plot(pca1$li[,1:2], col=color_ny, pch=19, cex=2, #$li is the principal components of the analysis;
     ylab=paste("PC2", paste("(", round(pca1$eig[2]/sum(pca1$eig)*100), " %)", sep="" )), #I don't understand the axis
     xlab=paste("PC1", paste("(", round(pca1$eig[1]/sum(pca1$eig)*100), " %)", sep="" )))
dev.off()
#
col <- Clusters$Cluster
rbPal <- colorRampPalette(c('red3','blue3'))
color_ny <- rbPal(20)[as.numeric(cut(col,breaks = 20))]
#
# For microsatellites -----------------------------------------------------

setwd('C:/Users/narjesy/Documents/Rworkspace/adegenet/')
ny_data <- read.table('./local/micro/local_microsatellites.txt', header = F, sep = '\t') #modified genepop file > no pop, no header
ny_data$V1 <- NULL
#as.character(ny_data) #Didn't change the plot
Clusters <- read.table('./local/micro/micro_clusters2.txt', header = T)
colnames(Clusters)
Clusters[1:2,]
#For colors:
col <- Clusters$Cluster
rbPal <- colorRampPalette(c('blue3','red3'))
color_ny <- rbPal(20)[as.numeric(cut(col,breaks = 20))]
#
#with adegenet

local_mic <- apply(ny_data, 2, function(x) as.character(x))
local_mic[local_mic=="0"] <- NA

library("adegenet")
library("ape")
local_mic_gen <- df2genind(local_mic, ploidy = 1, ncode = 3)

X <- scaleGen(local_mic_gen, NA.method="mean") #Replaces missing data with mean of allele frequencies
class(X)
dim(X)
#autoplot(prcomp(X), data= Clusters, colour = 'Clusters') #now it worked

##To get the axis values

pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3) #dudi.pco performs PCoA, it needs a distance matrix
#
plot(pca1$li[,1:2], pch=19, cex=2.5, col=color_ny, cex.lab=1.4,lyt=2, #$li is the principal components of the analysis;
     ylab=paste("PC2", paste("(", round(pca1$eig[2]/sum(pca1$eig)*100), " %)", sep="" )), #I don't understand the axis
     xlab=paste("PC1", paste("(", round(pca1$eig[1]/sum(pca1$eig)*100), " %)", sep="" )))

#For colors:
col <- Clusters$Clusters
rbPal <- colorRampPalette(c('blue3','red3'))
color_ny <- rbPal(20)[as.numeric(cut(col,breaks = 20))]
#
# For gcf data ------------------------------------------------------------
#This is from Petri, because the PC1 should be higher:
df_gcf <- read.table('./adegenet/gcf/74_p1r5maf05.genepop', header = F, sep = '\t') #modified genepop file > no pop, no header


df_gcf$V1 <- NULL #To get rid of the ind names. it seems like the software reads them as the alleles
df_gcf <- apply(df_gcf, 2, function(x) as.character(x))
df_gcf[5:5]
df_gcf[df_gcf=="0"] <- NA
gcf_clusters <- read.table('./adegenet/gcf/clusters.txt', header = T, sep = '\t') #modified genepop file > no pop, no header
colnames(gcf_clusters)
gcf_clusters[1:5,]

df_gcf_gen <- df2genind(df_gcf, ploidy = 1, ncode = 3)
#
col <- gcf_clusters$Cluster
rbPal <- colorRampPalette(c('blue3','green','red3'))
color_ny <- rbPal(20)[as.numeric(cut(col,breaks = 20))]
shapes <- gcf_clusters$Regions
#
#autoplot(prcomp(X), data= gcf_clusters, colour = 'Cluster') ##load ggfortify for this
##DAPC
grp <- find.clusters(df_gcf_gen, max.n.clust=10, n.pca=200) #NY: this was grp_NY in the code, I changed it to grp, bcs the next is grp...
3
str(grp)
is(grp)
length(grp[[3]])
dapc1 <- dapc(df_gcf_gen, grp$grp, n.pca = 3, n.da=2) #group by cluster
scatter(dapc1)
#PCA
X <- scaleGen(df_gcf_gen, NA.method="mean") #Replaces missing data with mean of allele frequencies
class(X)
dim(X)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3) #dudi.pco performs PCoA, it needs a distance matrix
#
s.label(pca1$li)
title("PCA of gcf magellanicum")
add.scatter.eig(pca1$eig[1:20], 3,1,2)
#
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50)) #Good!
tiff("./adegenet/gcf/gcf_adegenet_6-3.5_2.tiff", width = 6, height = 3.5, pointsize = 1/1000, units = 'in', res = 1000)
plot(pca1$li[,1:2], col=color_ny, pch = c(0, 2, 3, 1, 8)[as.numeric(shapes)], cex=2,lwd=3, #$li is the principal components of the analysis;
     ylab=paste("PC2", paste("(", round(pca1$eig[2]/sum(pca1$eig)*100), " %)", sep="" )), #I don't understand the axis
     xlab=paste("PC1", paste("(", round(pca1$eig[1]/sum(pca1$eig)*100), " %)", sep="" )))
dev.off()
#To add legend. I migt fix it later. It also needs a gradient...
legend("topleft",legend=paste(c("Asia","Eastern America","Eastern Canada","Europe", "Western Canada")),
       pch = c(0, 2, 3, 1, 8))
#

###To get the color gradient:

#

#Below didn't work
head(pca1)
plot(pca1$li[,1:2], col=grp$grp, pch=16)

ggplot(pca1, aes(x=li1, y=li2) + geom_point(size=3)) #Argggg, it doesn't work


# S. warnstorfii -----------------------------------------------------
getwd()
#PCA with genind made by vcfR2genind
df <- read.vcfR('M:/Project/R/Rworkspace_M/adegenet/warn/ref_fallax_m3_80/p1r5maf05/p1r5maf05_no_Ho.vcf', limit = 1e+07, nrows = -1, skip = 0, cols = NULL,
                convertNA = TRUE, verbose = TRUE)

df_genind <- vcfR2genind(df, sep = "[/]") 

#
warn_info <- read.table('M:/Project/R/Rworkspace_M/adegenet/warn/info_80.txt', header = T, sep = '\t') #
colnames(warn_info)
warn_info[1:5,1:5]
#
col <- warn_info$cluster_ref
rbPal <- colorRampPalette(c('red3','green','blue3'))
color_ny <- rbPal(20)[as.numeric(cut(col,breaks = 20))]
shapes <- warn_info$Region
levels(shapes)
#
#autoplot(prcomp(X), data= gcf_clusters, colour = 'Cluster') ##load ggfortify for this

#PCA
X <- scaleGen(df_genind, NA.method="mean") #Replaces missing data with mean of allele frequencies
class(X)
dim(X)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3) #dudi.pco performs PCoA, it needs a distance matrix
#
s.label(pca1$li)
title("PCA of S. warnstorfii")
add.scatter.eig(pca1$eig[1:20], 3,1,2)
#
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50)) #Good!
tiff("M:/Project/R/Rworkspace_M/adegenet/warn/ref_fallax_m3_80/p1r5maf05/p1r5maf05_no_Ho_2.tiff", width = 7, height = 3, pointsize = 1/1000, units = 'in', res = 1000)
par(mgp=c(1,0.4,0),mar=c(2,2,1,1)+0.1)
plot(pca1$li[,1], pca1$li[,2] ,  cex=2,lwd=3, col=color_ny , pch = c(2,1)[as.numeric(shapes)], #col=color_ny, # #$li is the principal components of the analysis; 
     cex.lab=0.8, cex.axis=0.8, 
     ylab=paste("PC2", paste("(", round(pca1$eig[2]/sum(pca1$eig)*100), " %)", sep="" )),
     xlab=paste("PC1", paste("(", round(pca1$eig[1]/sum(pca1$eig)*100), " %)", sep="" )))
#text(pca1$li[,1], pca1$li[,2], labels=warn_info$Country, cex= 0.4) # YAAAYYYY
legend("bottomleft",legend=paste(c('Arctic', 'Boreal', 'Temperate')),
     pch = c(2,1,3), cex=0.6)

dev.off()
#
# to plot PC 1 and 3, or PC2 and 3
tiff("M:/Project/R/Rworkspace_M/adegenet/warn/m3M0n1N0_72_blue/p1r5maf05_PC2_3.tiff", width = 4.5, height = 3.5, pointsize = 1/1000, units = 'in', res = 1000)
par(mgp=c(1.5,0.5,0),mar=c(3,3,2,2)+0.1)
plot(pca1$li[,2], pca1$li[,3] ,  cex=1,lwd=2, col= 'blue', pch = c(2,1,3)[as.numeric(shapes)], #col=color_ny, # #$li is the principal components of the analysis; 
     cex.lab=0.5, cex.axis=0.5, 
     ylab=paste("PC3", paste("(", round(pca1$eig[3]/sum(pca1$eig)*100), " %)", sep="" )), #I don't understand the axis
     xlab=paste("PC2", paste("(", round(pca1$eig[2]/sum(pca1$eig)*100), " %)", sep="" )))
legend("bottomright",legend=paste(c('Arctic', 'Boreal', 'Temperate')),
       pch = c(2,1,3), cex=0.4)
dev.off()
#
##DAPC
grp <- find.clusters(df_genind, max.n.clust=10, n.pca=200) #NY: this was grp_NY in the code, I changed it to grp, bcs the next is grp...
2
str(grp)
is(grp)
length(grp[[3]])
dapc1 <- dapc(df_genind, grp$grp, n.pca = 3, n.da=2) #group by cluster
scatter(dapc1)

#


#

# 



