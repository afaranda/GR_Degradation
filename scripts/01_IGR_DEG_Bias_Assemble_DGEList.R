##############################################################################
#                                                                            #
# File: 01_IGR_DEG_Bias_Assemble_DGEList.R                                   #
# Puropose: Prepare gene level count data from SRA Project SRP041955 for     #
#           analysis and create a corresponding synapse project. This script #
#           organizes gene level counts, gene and sample metadata into a     #
#           DGEList object that can be used to run downstream analyses       #
#                                                                            #
# Created: July 9, 2021                                                      #
# Author: Adam Faranda                                                       #
#                                                                            #
##############################################################################

############ Load Libraries and Source Dependencies ##########################

library('synapser')
library('dplyr')
library('edgeR')

wd <- "/Users/adam/Documents/31Aug2021_Degradation_Bias_Assesment_in_SRP041955/"
setwd(wd)
local_data_dir <- paste0(wd,"/data")           ## Local data directory
local_res_dir <-  paste0(wd,"/results")        ## Local results directory
########################## Setup Synapse Project #############################
synapser::synLogin()
syn_project <- Project("Length Bias in Deliberately Degraded Samples")
syn_project <- synStore(syn_project)

# Create folder for input data
syn_count_dir <- Folder(
  "data",
  parent=syn_project
)
syn_count_dir <- synStore(syn_count_dir)

# Create folder to store code
syn_code_dir <- Folder(
  "code",
  parent=syn_project
)
syn_code_dir <- synStore(syn_code_dir)

# Create folder to results 
syn_res_dir <- Folder(
  "results",
  parent=syn_project
)
syn_res_dir <- synStore(syn_res_dir)

# Setup sample covariate table and push count files to Synapse
count_files <- list.files(local_data_dir)
ft <- data.frame()
used_list <- list()  # Start a list of data files used by this script
for(f in count_files) {
  if(grepl("_rf_GeneCount", f)){
    ft <- rbind(
      ft, 
      data.frame(
        sample=gsub("_rf_GeneCount.txt", "", f),
        files=f
      )
    )
    synCountFile <- File(
      paste(local_data_dir, f,sep="/"), 
      parent=syn_count_dir
    )
    synCountFile <- synStore(synCountFile)
    used_list <- append(used_list, synCountFile)
  }
}

# Add covariates to sample table
ft <- ft%>%
  mutate(
    subject=gsub(
        "_[0123456789]*","",
        gsub("SAMN[0123456789]{8}_","", sample)
      ),
    interval=factor(
      paste0(
        gsub("SAMN[0123456789]{8}_[1234]_","", sample),"H"
      ), levels=c("0H", "12H", "24H", "48H", "84H")
    )
  ) %>%
  mutate(
    label=paste0("Subj_",subject,"_", interval)
  ) %>%
  tibble::column_to_rownames("label")%>%
  mutate(
    label=paste0("Subj_",subject,"_", interval)
  ) 

# Read count files into a dgelist object and store
master <- readDGE(
  files=ft,
  header=F, path=local_data_dir
)

# Strip QC Rows from Count Matrix
master <- master[!grepl("__", row.names(master$counts)),,keep.lib.sizes=F]

## add gene annotations to master DGEList, push annotation file to synapse
fn<-'data/Gene_Annotations.csv'
lt<-read.csv(fn, row.names = 1)
row.names(lt) <- lt$gene_id
master$genes <- lt[row.names(master$counts),]
synAnnot <- File(fn,parent=syn_count_dir)
synAnnot <- synStore(
  synAnnot
)


# Add this script to the code dir and update provenance
synapse_push <- File(
  path="scripts/01_IGR_DEG_Bias_Assemble_DGEList.R",
  parent=syn_code_dir
)

synapse_push <- synStore(
  synapse_push,
  used=append(used_list, synAnnot)
)

save(master, file="data/master_dgelist.Rdata")
syn_dge <- File(
  path="data/master_dgelist.Rdata",
  parent=syn_count_dir
)
syn_dge <- synStore(
  syn_dge
)

syn_act <- Activity(
  name="upload_dgelist",
  description="retrieved counts from synapse and assembled a DGEList"
)

syn_act$executed(synapse_push)

synSetProvenance(
  syn_dge,
  syn_act
)


