library(dplyr)
library(data.table)
library(rMVP)
MVP.Data(fileVCF = 'SbDiv3458chr.vcf',
         filePhe = 'pheno.csv',
         sep.phe = ',',
         fileKin = TRUE,
         filePC = TRUE,
         out = 'mvp.vcf')MVP.Data.VCF2MVP("SbDiv_imputed_Filtermaf.recode.vcf", out='mvp')MVP.Data.PC("MVP.pc.desc", out='MVP', sep='\t')
MVP.Data.PC(TRUE, mvp_prefix='MVP', perc=1, pcs.keep=5)
pca <- attach.big.matrix("mvp.vcf.pc.desc")
pca2 <- attach.big.matrix("mvp.vcf.pc.desc")[, 1:3]
MVP.PCAplot(PCA=pca, Ncluster=3, class=NULL, col=c("red", "green", "yellow"), file.type="jpg")
write.csv(pca2, "pcaanalysis.csv", row.names = FALSE)
