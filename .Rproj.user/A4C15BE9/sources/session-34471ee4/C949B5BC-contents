##############################################################################
# File: 04_DEG_Analysis.R                                                    #
# Purpose:                                                                   #
#       Examine the correnation between log2 fold change and gene coding     #
#       length in data from Irene Gallego Romero's 2014 paper investigating  #
#       differential transcript stability in human blood samples maintained  #
#       at room temperature for intervals between 0 and 84 hours after       #
#       collection and prior to library preparation. Blood samples drawn     #
#       from four human subjects (subject 1-4).                              #
# Author: Adam Faranda                                                       #
##############################################################################

############################ Setup Environment ###############################

# CHANGE THIS DIRECTORY
wd <- "/Users/adam/Documents/31Aug2021_Degradation_Bias_Assesment_in_SRP041955/"
setwd(wd)
local_data_dir <- paste0(wd,"/data")           ## Local data directory
local_res_dir <-  paste0(wd,"/results")        ## Local results directory

library(dplyr)
library(cluster)
library(reshape2)
library(pheatmap)
library(sva)
library(synapser)
source('scripts/02_Wrap_DEG_Functions.R')
source('scripts/03_Overlap_Comparison_Functions.R')

synapser::synLogin()
syn_project <-synFindEntityId("Length Bias in Deliberately Degraded Samples")

# Get Synapse ID of folder for input data
syn_count_dir <- synFindEntityId(
  "data",
  parent=syn_project
)

# Get Synapse ID of folder for code
syn_code_dir <- synFindEntityId(
  "code",
  parent=syn_project
)

# Get Synapse ID of folder for results 
syn_res_dir <- synFindEntityId(
  "results",
  parent=syn_project
)

########################### Fetch Master DGE List ############################

fn <- "master_dgelist.Rdata"
syn_dge <- synFindEntityId(name=fn, parent=syn_count_dir)
if(file.exists(paste0(local_data_dir,"/",fn))){
  print("Loading Local DGEList")
  load(paste0(local_data_dir,"/",fn))
} else if(!is.null(syn_dge)){
  print("Fetching DGEList from Synapse")
  synGet(
    syn_dge,
    downloadLocation = local_data_dir
  )
} else {
  print("Assembling DGEList from Counts")
  source("scripts/01_IGR_DEG_Bias_Assemble_DGEList.R")
}

######################### Analyze Wildtype Samples ###########################
# DBI Analysis ####
dge<-master[,master$samples$subject != 2]   # Select Samples
dge$samples$group <- dge$samples$interval



# Define Experimental Design
design<-model.matrix(~0+interval, dge$samples)     
# colnames(design)<-gsub(
#   'interval', '', 
#   colnames(design)
# )

# Define Contrasts
cntmat<-makeContrasts(                         
  Degrade_12H = `interval12H` - `interval0H`,
  Degrade_24H = `interval24H` - `interval0H`,
  Degrade_48H = `interval48H` - `interval0H`,
  Degrade_84H = `interval84H` - `interval0H`,
  levels = design
)

# Normalize DGEList and fit model
obj <- process_edgeR_ByDesign(
  dge,
  design=design
)

# Create Diagnostic Figures
diagnostic_plots(
  dge=obj$dge, prefix = "DBI_Wildtype",
  respath = local_res_dir
)

## Iterate over contrasts and prepare summaries / DEG tables
dge$samples$group <- paste0("interval",dge$samples$interval)
obj <- process_edgeR_ByDesign(
  dge,
  design=design
)

df <- data.frame()
deg <- data.frame()
res <- iterate_edgeR_pairwise_contrasts(
  obj[[1]], obj[[2]], cntmat, df=df, design=design,
  deg=deg, respath = local_res_dir,
  prefix = "IGR_Deg_Len_Bias"
)

######## Plot Bias for Statistically and Biologiclly significant DEG #########
length_bias_plot <- function(deg, fn){
  Contrast <- paste0(
    unique(deg$Group_2,), "v",
    unique(deg$Group_1)
  )
  
  mn<-paste("Length Bias in",Contrast)
  yl <-paste(                                    # Construct y axis label
    "Log2 Fold Change in",  unique(deg$Group_2,),          
    "vs", unique(deg$Group_1,)
  )
  
  png(fn, width=6, height = 5, units="in", res=1200)
  plot(
    x=log(deg$eu_length,2), y=deg$logFC, main=mn,
    xlab = "Log2 Gene Length (exon-union)",
    ylab = yl,
    col = ifelse(
      deg[,"Avg1"] == 0,
      "red",
      ifelse(deg[,"Avg2"] == 0, "red", "black")
    )
  )
  
  # Add Correlation Coefficients to tests
  ct=cor.test(deg$logFC, log(deg$eu_length,2), method="spearman")
  rho=paste("rho:", round(ct$estimate,3))
  sig=paste("p value:", signif(ct$p.value,3), sep="")
  abline(lm(deg$logFC~log(deg$eu_length,2)), col="red", lwd=2)
  text(7,max(deg$logFC)-1,rho)
  text(7,max(deg$logFC)-3,sig)
  dev.off()
}

gc_bias_plot <- function(deg, fn){
  Contrast <- paste0(
    unique(deg$Group_2,), "v",
    unique(deg$Group_1)
  )
  
  mn<-paste("GC Bias in",Contrast)
  yl <-paste(                                    # Construct y axis label
    "Log2 Fold Change in",  unique(deg$Group_2,),          
    "vs", unique(deg$Group_1,)
  )
  
  png(fn, width=6, height = 5, units="in", res=1200)
  plot(
    x=deg$eu_gc, y=deg$logFC, main=mn,
    xlab = "Fractional GC content (exon-union)",
    ylab = yl,
    col = ifelse(
      deg[,"Avg1"] == 0,
      "red",
      ifelse(deg[,"Avg2"] == 0, "red", "black")
    )
  )
  
  # Add Correlation Coefficients to tests
  ct=cor.test(deg$logFC, deg$eu_gc, method="spearman")
  rho=paste("rho:", round(ct$estimate,3))
  sig=paste("p value:", signif(ct$p.value,3), sep="")
  abline(lm(deg$logFC~deg$pl_gc), col="red", lwd=2)
  text(0.35,max(deg$logFC)-1,rho)
  text(0.35,max(deg$logFC)-3,sig)
  dev.off()
}


deg <- res[[2]] %>%
  mutate(
    Group_1=gsub("interval", "", Group_1),
    Group_2=gsub("interval", "",Group_2)
  )

for(c in unique(deg$Group_2)){
  fn <- paste0(
    local_res_dir,
    "/IGR_Length_Bias_",c,"_vs_0H_Stat.png"
  )
  
  length_bias_plot(deg %>% filter(Group_2 == c & Test=="ExactTest"),fn)
  
  fn <- paste0(
    local_res_dir,
    "/IGR_Length_Bias_",c,"_vs_0H_Bio.png"
  )
  length_bias_plot(
    deg %>% filter(
      Group_2 == c &
        Test=="ExactTest" &
        abs(Avg1 - Avg2) > 2 &
        (Avg1 > 2 | Avg2 > 2)
    ), 
    fn
  )
}

res[[2]] %>% group_by(Test, Group_1, Group_2) %>%
  summarize(Bias = cor(logFC, log2(eu_length), method="spearman"))

res[[1]] %>% filter(criteria == "Statistically Significant") %>%
  arrange(test)

res[[2]] %>% group_by(Test, Group_1, Group_2) %>%
  filter(abs(Avg1 - Avg2) > 2 & (Avg1 > 2 | Avg2 > 2))%>%
  summarize(Bias = cor(logFC, log2(eu_length), method="spearman"))


res[[1]] %>% filter(criteria == "Biologically Significant") %>%
  arrange(test)



###################### Push Analysis Results to Synapse ######################
# Main Analysis Script
syn_main_script <- File(
  path="transcriptomic_analysis_scripts/LIRTS_DEG_Analysis.R",
  parent=syn_code_dir
)
syn_main_script <- synStore(
  syn_main_script
)

# Wrappers for DEG Analysis Functions
syn_wrap_script <- File(
  path="transcriptomic_analysis_scripts/LIRTS_Wrap_DEG_Functions.R",
  parent=syn_code_dir
)
syn_wrap_script <- synStore(
  syn_wrap_script
)

# Helper script for deg summaries
syn_ovl_script <- File(
  path="transcriptomic_analysis_scripts/Overlap_Comparison_Functions.R",
  parent=syn_code_dir
)
syn_ovl_script <- synStore(
  syn_ovl_script
)

## Construct activity to set provenance on
syn_act <- Activity(
  name="DEG_Analysis",
  description=paste(
    "Analyze various contrasts in the LIRTS",
    "data set for differential expression"
  ),
  executed=c(
    syn_main_script,
    syn_wrap_script,
    syn_ovl_script
  ),
  used=syn_dge
)

# Add a LIRTS_DEG_Tables Directory to store individual DEG Tables
syn_fpkm_dir <- Folder(name="LIRTS_Expression_Matrices", parent=syn_project)
syn_fpkm_dir <- synStore(
  syn_fpkm_dir
)

for(f in list.files("LIRTS_DEG_Analysis_results/", pattern="FPKM_Matrix")){
  syn_fpkm <- File(
    path=paste0(
      "LIRTS_DEG_Analysis_results/", f
    ),
    parent=syn_fpkm_dir
  )
  synStore(syn_fpkm, activity = syn_act)
}

# Add Feature Selection tables to the Genes folder
for(
  f in list.files(
    "LIRTS_DEG_Analysis_results/", 
    pattern="Feature_Selection_Table")
){
  syn_fst <- File(
    path=paste0(
      "LIRTS_DEG_Analysis_results/", f
    ),
    parent=syn_gene_meta
  )
  synStore(syn_fst, activity = syn_act)
}

# Add a LIRTS_DEG_Tables Directory to store individual DEG Tables
syn_deg_dir <- Folder(name="LIRTS_DEG_Tables", parent=syn_project)
syn_deg_dir <- synStore(
  syn_deg_dir
)

for(f in list.files("LIRTS_DEG_Analysis_results/", pattern="Test_DEG")){
  syn_deg <- File(
    path=paste0(
      "LIRTS_DEG_Analysis_results/", f
    ),
    parent=syn_deg_dir
  )
  synStore(syn_deg, activity = syn_act)
}

# Add a LIRTS_DEG_Tables Directory to store Princapal Components and
# BCV Plots (Diagnostic Plots)
syn_fig_dir <- Folder(name="LIRTS_DEG_PCA_BCV_Plots", parent=syn_project)
syn_fig_dir <- synStore(
  syn_fig_dir
)

for(f in list.files("LIRTS_DEG_Analysis_results/", pattern="PCA_Plot")){
  syn_pca <- File(
    path=paste0(
      "LIRTS_DEG_Analysis_results/", f
    ),
    parent=syn_fig_dir
  )
  synStore(syn_pca, activity = syn_act)
}

for(f in list.files("LIRTS_DEG_Analysis_results/", pattern="BCV_Plot")){
  syn_pca <- File(
    path=paste0(
      "LIRTS_DEG_Analysis_results/", f
    ),
    parent=syn_fig_dir
  )
  synStore(syn_pca, activity = syn_act)
}

## Create master DEG table in Synapse (if none exists)
syn_deg_table <- synFindEntityId(
  "LIRTS_DEG_Master_Table", 
  parent=syn_project
)

if(is.null(syn_deg_table)){
  syn_cols <- list(                                 
    Column(name="gene_id", columnType='STRING', maximumSize =100), 
    Column(name='logFC', columnType='DOUBLE'),
    Column(name='logCPM', columnType='DOUBLE'),
    Column(name='PValue', columnType='DOUBLE'),
    Column(name='FDR', columnType='DOUBLE'),
    Column(name='Avg1', columnType='DOUBLE'),
    Column(name='Avg2', columnType='DOUBLE'),
    Column(name='Group_1', columnType='STRING', maximumSize=50),
    Column(name='Group_2', columnType='STRING', maximumSize=50),
    Column(name='Test', columnType='STRING', maximumSize=100),
    Column(name='Samples', columnType='STRING', maximumSize=100)
  )
  syn_schema <- Schema(
    name="LIRTS_DEG_Master_Table",
    columns=syn_cols,
    parent=syn_project
  )
  syn_table <- Table(
    schema = syn_schema,
    values=res[[2]]
  )
  
  syn_table <- synStore(syn_table, activity = syn_act)
}