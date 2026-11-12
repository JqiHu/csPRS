# Assign SNPs in GWAS to genes in each cells
setwd('/path/to/csPRS/score/weights_by_cell_type')
library(stringr)
library(data.table)
library(purrr)

# GWAS summary statistics with position and weights from SDPR
### cleaned GWAS
#gwas <- readRDS('/path/to/resources/prs/sumstats/cleaned_sumstats/cad_cardiogramc4a_sum_rm_dup_overlapUKB.RDS')
### SDPR weights
#weight <- fread('/path/to/resources/prs/beta/SDPR_cad_cardiogramc4a.txt',data.table=F)

#dat <- merge(weight,gwas,by.x='SNP',by.y='V2')

## hg19 to hg38
#bim <- subset(dat,select=c('CHR','POS','SNP'))
#colnames(bim) <- c('chr','pos','SNP')
#library(bigsnpr)
#liftover <- '/path/to/software/liftOver'
#i=T
#while(i){
#  try_err <- try(bim38 <- snp_modifyBuild(bim, liftover, 
#                                          from='hg19', to='hg38'))
#  if(class(try_err)[1]!='try-error'){i=F}
#}
## 513 missing
#colnames(bim38)[2] <- 'pos_hg38'
#dat2 <- merge(bim38,dat,by.x=c('chr','SNP'),by.y=c('CHR','SNP'))

#write.table(dat2,gzfile('../CAD_clean_hg38.txt.gz'),
#	    row.names=F,sep='\t',quote=F)

#dat2 <- fread('../CAD_clean_hg38.txt.gz',data.table=F)

#dat2 <- dat2[is.na(dat2$pos_hg38)==F,]
# Function to extract SNPs in certain region
#getSNPs <- function(chr,start,stop,gene,gene_name){
#  out <- dat2[dat2$chr==chr &
#	     dat2$pos_hg38>=(start-10000) &
#	     dat2$pos_hg38<=(stop+10000),]
#  res <- out
#  if(nrow(res)>0){
#    res$gene <- gene
#    res$gene_name <- gene_name

#    print(nrow(res))
#    return(res)
#  }
#}
### Extract SNPs for each gene
#genes <- fread('../../gene_list/gene_info.csv',data.table=F)
#colnames(genes)[c(1,2,4:6)] <- c('gene','gene_name','chr','start','stop')
#genes2 <- genes[,c(1,2,4:6)]

#snp_to_gene <- pmap(genes2, function(chr, start, stop, gene, gene_name) {
#                          getSNPs(chr, start, stop, gene, gene_name)
#	    })
#names(snp_to_gene) <- genes2$gene

#saveRDS(snp_to_gene,'snp_to_gene.rds')
snp_to_gene <- readRDS('snp_to_gene.rds')

# Get results cell by cell
#path1 <- '../../gene_list/genes_by_cell_type/'
#files <- list.files(path1)
#all <- c()
#for(f in files){
#  tmp <- fread(paste0(path1,f),data.table=F) # hg38!!!!!!!!
#  res <- snp_to_gene[tmp$ensembl_gene_id]
  
#  add <- do.call(rbind,res)
#  add <- add[duplicated(add$SNP)==F,]
#  print(dim(add))
#  cell <- str_remove_all(f,'_genes.txt')
#  add$cell <- cell

#  write.table(add,paste0(cell,'_weights.txt'),
#	      row.names=F,sep='\t',quote=F)
  
#  all <- rbind(all,add)
#  print(paste0('Finished for ',f))
#}
### save all the cells
#write.table(all,gzfile('SDPR_CAD_cell_type.txt.gz'),
#	    row.names=F,sep='\t',quote=F)

### Count number of genes in each cell
#count1 <- aggregate(all$gene_name,list(all$cell),
#		    function(x){length(unique(x))})
### Count number of SNPs in each cell
#count2 <- aggregate(all$SNP,list(all$cell),
#		    function(x){length(unique(x))})
#count <- merge(count1,count2,by='Group.1')
#colnames(count) <- c('cellType','num_gene','num_SNP')
#write.csv(count,'../description/count_snp_cell_type.csv',row.names=F)


# Get results cell by cell
setwd('../weights_by_cell_state')
path2 <- '../../gene_list/genes_by_cell_state/'
files <- list.files(path2)
all <- list()
for(f in files){
  tmp <- fread(paste0(path2,f),data.table=F) # hg38!!!!!!!!
  res <- snp_to_gene[tmp$ensembl_gene_id]

  add <- do.call(rbind,res)
  add <- add[duplicated(add$SNP)==F,]
  print(dim(add))
  cell <- str_remove_all(f,'_genes.txt')
  add$cell <- cell

  write.table(add,paste0(cell,'_weights.txt'),
              row.names=F,sep='\t',quote=F)

  all[[cell]] <- add
  print(paste0('Finished for ',f))
}
### save all the cells
saveRDS(all,'SDPR_CAD_cell_state.rds')

### Count number of genes in each cell
count1 <- lapply(all,function(x){length(unique(x$gene_name))})
### Count number of SNPs in each cell
count2 <- lapply(all,function(x){nrow(x)})

count <- data.frame(cellType=names(all),
		    num_gene=unlist(count1),
		    num_SNP=unlist(count2))
print(count)

write.csv(count,'../description/count_snp_cell_state.csv',row.names=F)


