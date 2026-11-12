# organize sLDSC results
setwd('/path/to/csPRS/heritability_par/res')
library(data.table)
suffix <- 'merge251127'

# Read results
res1 <- fread(paste0('CAD_',suffix,'.results'),data.table=F)
### keep only cells
res1 <- res1[1:(grep('baseL2_1',res1$Category)-1),]
#res2 <- fread(paste0('CAD_cell_state_',suffix,'.results'),data.table=F)
#res2 <- res2[1:(grep('baseL2_1',res2$Category)-1),]

res <- res1
#res <- rbind(res1,res2)

### p-value threshold
thresh <- 0.05/nrow(res)
res.sig <- res[res$Enrichment_p<thresh,]
res.sig <- res.sig[is.na(res.sig$Category)==F,] # nothing significant

write.csv(res,paste0('CAD_cell_all_',suffix,'.csv'),row.names=F)
