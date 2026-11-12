# Create weights for top 1000 genes
setwd('/path/to/csPRS/score/weights/')
library(data.table)
library(stringr)

## Function to generate outputs by types
getOut <- function(cell,path_out,path_weight,path_gene){
  # Read weights
  weight <- fread(paste0(path_weight,cell,'_weights.txt'),data.table=F)

  # Read genes
  genes <- readRDS(paste0(path_gene,'top10_gene.rds'))
  names(genes) <- str_replace_all(names(genes),'_','-')
  genes_cell <- genes[[cell]]

  out <- weight[weight$gene_name %in% genes_cell,]
  print(nrow(out))

  write.table(out,paste0(path_out,cell,'.txt'),
              row.names=F,sep='\t',quote=F)
}

# ---------------
# For cell types
## all weight files
files <- list.files('weights_by_cell_type','weights.txt')
types <- str_remove_all(files,'_weights.txt')
## path of output
path_out='weights_specificity_top10pct_cell_type/'
## path to overall weight
path_weight='weights_by_cell_type/'
## path to genes
path_gene='/path/to/csPRS/gene_list/genes_specitificy_top10/cell_type/'

### Run for all types
lapply(types,function(x){getOut(x,path_out,path_weight,path_gene)})

# ---------------
# For cell states
## all weight files
files <- list.files('weights_by_cell_state','.txt')
types <- str_remove_all(files,'_weights.txt')
## path of output
path_out='weights_specificity_top10pct_cell_state/'
## path to overall weight
path_weight='weights_by_cell_state/'
## path to genes
path_gene='/path/to/csPRS/gene_list/genes_specitificy_top10/cell_state/'

### Run for all types
lapply(types,function(x){getOut(x,path_out,path_weight,path_gene)})

