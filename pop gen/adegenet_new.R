getwd()
setwd('/Users/narcis/Narcis/Project/Hottonia/palustris/ddRAD/Mallorca/')

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

df <- read.vcfR('./adegenet/pop_r.7_p2/populations.snps.vcf', limit =20000000, nrows = -1, skip = 0, cols = NULL, convertNA = TRUE, verbose = TRUE)

dp <- extract.gt(df, element = "DP", as.numeric=TRUE)
hist(dp)
?hist

# to get percentage of missing data
myMiss <- apply(dp, MARGIN = 2, function(x){ sum(is.na(x)) })
myMiss <- myMiss/nrow(df)

# to plot it
library(RColorBrewer)
palette(brewer.pal(n=12, name = 'Set3'))

par(mar = c(12,4,4,2))
barplot(myMiss, las = 2, col = 1:12)
title(ylab = "Missingness (%)")
par(mar = c(5,4,4,2))
#
# make genind object
df_genind <- vcfR2genind(df, sep = "[/]") 

#NB: another way to read data in is import2genind() from adegenet. it tries to detect the file format and finds the appropriate import function
# for import2genind(), structure should be with .str or .stru, and genepop with .gen extention
# you can also use read.structure or read.genepop

#PCA
X <- scaleGen(df_genind, NA.method="mean") #Replaces missing data with mean of allele frequencies
class(X)
dim(X)
pcoa1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3) #dudi.pco performs PCoA, it needs a distance matrix
#
s.label(pcoa1$li)
title("PCA of H. palustris")
add.scatter.eig(pcoa1$eig[1:9], 3,1,2)
#
popmap <- read.table('./pop_gen/popmap_2', header = T)
names(popmap)
pops <- popmap$populations
#

barplot(pcoa1$eig[1:10],main="PCA eigenvalues", col=heat.colors(50)) #Good!
tiff("./adegenet/barplot.tiff", width = 7.5, height = 3.5, pointsize = 1/1000, units = 'in', res = 1000)
par(mgp=c(3,0.8,0),mar=c(4,4,1,1)+0.1)
plot(pcoa1$li[,1], pcoa1$li[,2], cex=2, lwd=3, pch = 1, #col=color_ny, # #$li is the principal components of the analysis; 
     cex.lab=1, cex.axis=1, 
     ylab=paste("PC2", paste("(", round(pcoa1$eig[2]/sum(pcoa1$eig)*100), " %)", sep="" )),
     xlab=paste("PC1", paste("(", round(pcoa1$eig[1]/sum(pcoa1$eig)*100), " %)", sep="" )),
     text(pcoa1$li[,1], pcoa1$li[,2], labels=popmap$ind, cex= 0.7)) 
#legend("bottomright",legend=paste(c('France', 'Germany', 'Switzerland')),
       #pch = c(2,1,3), cex=0.6)

dev.off()

#
#
#
df_genlight <- vcfR2genlight(df, n.cores = 2) 
??vcfR2genlight
??glPca
??glDotProd
pca1 <- glPca(df_genlight)

# everything related to genlight class is following http://adegenet.r-forge.r-project.org/files/tutorial-genomics.pdf
# df_genlight <- vcfR2genlight(df, parallel = TRUE) 
#
pca1 <- glPca(df_genlight)
pca1
scatter(pca1, posi="topleft") 
title("PCA of H. palustris\n axes 1-2")
# to add color
myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4) 
abline(h=0,v=0, col="grey") 
add.scatter.eig(pca1$eig[1:9],2,1,2, posi="topright", inset=.17, ratio=.16)

# to confirm the pca with NJ tree:
library(ape) 
tre <- nj(dist(as.matrix(df_genind))) 
tre

plot(tre, typ="fan", cex=0.7) 
title("NJ tree of H. palustris")
?plot
# tree with color
plot(tre, typ="fan", show.tip=FALSE) 
tiplabels(pch=20, col=myCol, cex=4) 
title("NJ tree of H. palustris")

# DAPC
#
#
Hpal_info <- read.table('./pop_gen/popmap_2', header = T)
colnames(Hpal_info)
#
col <- Hpal_info$cluster
rbPal <- colorRampPalette(c('red3','green', 'blue3')) # 'blue3','green','red3'
color_ny <- rbPal(20)[as.numeric(cut(col,breaks = 20))]
shapes <- Hpal_info$populations
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
title("PCA of H.palustris")
add.scatter.eig(pca1$eig[1:9], 3,1,2)
#
#barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50)) #Good!
tiff("./adegenet/Rplot_pca2.png", width = 7.5, height = 3.5, pointsize = 1/1000, units = 'in', res = 1000)
par(mgp=c(3,0.8,0),mar=c(4,4,1,1)+0.1)

plot(pca1$li[,1], pca1$li[,2] ,  cex=2,lwd=3, col=color_ny , pch = c(2,1,3)[as.numeric(shapes)], #col=color_ny, # #$li is the principal components of the analysis; 
     cex.lab=1, cex.axis=1, 
     ylab=paste("PC2", paste("(", round(pca1$eig[2]/sum(pca1$eig)*100), " %)", sep="" )),
     xlab=paste("PC1", paste("(", round(pca1$eig[1]/sum(pca1$eig)*100), " %)", sep="" )))
#text(pca1$li[,1], pca1$li[,2], labels=warn_info$Country, cex= 0.4) # YAAAYYYY
legend("bottomright",legend=paste(c('France', 'Germany', 'Switzerland')),
       pch = c(2,1,3), cex=0.6)

dev.off()
#
##DAPC
grp <- find.clusters(df_genind, max.n.clust=5, n.pca=200)
3
str(grp)
is(grp)
length(grp[[3]])
dapc1 <- dapc(df_genind, grp$grp, n.pca = 3, n.da=2) #group by cluster
scatter(dapc1)
#
#


