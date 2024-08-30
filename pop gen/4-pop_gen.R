# pop gen analysis

library(pegas)

??pegas

library(help = pegas)

# my_file <- read.vcf('/Users/narcis/Narcis/Project/Hottonia/palustris/ddRAD/analysis/Irina/populations.snps.vcf', from = 1, to = 10000, which.loci = NULL, quiet = FALSE)
df <- read.vcf('./pop_gen/pop_r.7_p2/Hpal_pop_r.7_p2.vcf', from = 1, to = 10000, which.loci = NULL, quiet = FALSE)
#
hw <- hw.test(df, B = 1000)
dim(hw)

#?write.table
write.table(hw, file = './hw.txt', sep = "\t")

#
popmap <- read.table('./pop_gen/popmap_2', header = T)
names(popmap)
pops <- popmap$populations


fst <- Fst(df, pop = pops, quiet = TRUE, na.alleles = "")
dim(fst)
head(fst)

hist(fst[,2])

write.table(fst, file = './fst.txt', sep = "\t")

#
het <- heterozygosity(df, variance = FALSE)
hist(het[,1])

write.table(het, file = './het.txt', sep = "\t")
#

## with package hierfstat. check here: https://cran.r-project.org/web/packages/hierfstat/vignettes/hierfstat.html

# install.packages('hierfstat')
library(hierfstat)
??hierfstat

df2 <- read.vcfR('./adegenet/pop_r.7_p2/populations.snps.vcf', limit =20000000, nrows = -1, skip = 0, cols = NULL, convertNA = TRUE, verbose = TRUE)
df_genind <- vcfR2genind(df2, sep = "[/]") 
df_hier <- genind2hierfstat(df_genind, pop=pops)

#df_dos <- fstat2dos(df_hier,diploid=TRUE)

basic_stats <- basic.stats(df_hier,diploid=TRUE,digits=4)
Fis <- basic_stats$Fis
write.table(basic_stats$Fis, file = "./basic_stats_Fis.txt", sep = "\t")
#
head(basic_stats$perloc)
write.table(basic_stats$perloc, file = "./basic_stats_perloc.txt", sep = "\t")

#boxplot
boxplot(basic_stats$perloc[,1:10])

#
basic_stats_overal <- basic_stats$overall
#
pairwise.neifst(df_hier,diploid=TRUE)

#
allele_count <- allele.count(df_hier,diploid=TRUE)
head(allele_count$X20.53)
allele_count$X20.53

allelic_richness <- allelic.richness(df_hier,min.n=NULL,diploid=TRUE)
head(allelic_richness$Ar)


g.stats(df_hier,diploid=TRUE)

basic_stats <- basic.stats(df_hier,diploid=TRUE,digits=4)
Fis <- basic_stats$Fis
write.table(basic_stats$Fis, file = "./basic_stats_Fis.txt", sep = "\t")
head(basic_stats$perloc)

pairwise.neifst(df_hier,diploid=TRUE)

x <- indpca(df_hier,ind.labels=NULL,scale=FALSE)
plot(x,cex=0.7, pch = 1)

#
distance <- genet.dist(df_hier,diploid=TRUE,method="Dch")
#
# pp.fst(dat=df_hier,diploid=TRUE) # Fst per pair

# test.between(data = df_hier, test.lev = pops, rand.unit = 2, nperm = 1000) # didn't work
?test.between.within()

?test.between

?diploid

