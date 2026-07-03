library(data.table)
library(dplyr)
library(tidyr)
library(parallel)
library(ggplot2)
library(GenomicRanges)
library(reticulate)

#Needed: ages_data file containing a SNP column, created like this:
#ages_data <- fread("Selection_Summary_Statistics_01OCT2025.tsv.gz")
#ages_data$SNP <- paste(ages_data$CHROM, ages_data$POS, ages_data$REF, ages_data$ALT, sep = ":")
#Call it ages_data

#Write names of all traits
traits <- c("urate", "LDL-C")

#Write file names of all GWAS summary stats
gwas_files <- c("30880_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz", "30780_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz")

for (i in seq_along(traits)) {
  
  trait <- traits[i]  
  
  cat("annotating", trait, "clumps", "\n")
  
  #Read in gwas summary stats file for trait 1  
  gwas <- fread(gwas_files[i])
  
  #Create separate columns for chromosome, position, reference and alternative alleles
  gwas <- gwas %>% separate(variant, into = c("chromosome", "position", "reference", "alternative"), sep = ":")
  #Create a SNP column again
  gwas$SNP <- paste(gwas$chromosome, gwas$position, gwas$reference, gwas$alternative, sep = ":") 
  
  #Make sure that str of both lists/dfs is the same
  gwas$chromosome <- as.integer(gwas$chromosome)
  gwas$position <- as.integer(gwas$position)
  gwas$SNP <- as.character(gwas$SNP)
  gwas$reference <- as.character(gwas$reference)
  gwas$alternative <- as.character(gwas$alternative)
  
  #Make a list of all snps in df_all_chr
  load(paste0(gwas_files[i],".rdat"))
  snps_in_df_all_chr <- unique(c(df_all_chr$snp1, df_all_chr$snp2))
  
  #Keep only the rows in LD_positions for which the snps are in df_all_chr
  load(paste0(gwas_files[i],"_clumps.rdat"))
  LD_positions_new <- lapply(LD_positions, function(df) {
    df[df$rsid %in% snps_in_df_all_chr, ]
  })
  
  #Filter gwas data for significant variants
  gwas_sign <- subset(gwas, pval < 5e-8)
  
  #Keep only the snps that are significant in the GWAS (by merging)
  LD_positions_sign <- lapply(LD_positions_new, function(df) {
    merge(df, gwas_sign[, c("chromosome", "position", "reference", "alternative")],
          by.x = c("chromosome", "position", "allele1", "allele2"),
          by.y = c("chromosome", "position", "reference", "alternative"),
          all = FALSE,
          sort = FALSE)
  })
  
  #Make sure that the column names are the same for all dataframes
  LD_positions_sign <- lapply(LD_positions_sign, function(df) {
    names(df)[names(df) == "allele1"] <- "reference"
    df
  })
  
  LD_positions_sign <- lapply(LD_positions_sign, function(df) {
    names(df)[names(df) == "allele2"] <- "alternative"
    df
  })
  
  #Create a SNP column for each df in the list
  LD_positions_sign <- lapply(LD_positions_sign, function(df) {
    df$SNP <- paste(df$chromosome, df$position, df$reference, df$alternative, sep = ":")
    df
  })
  
  #Add the necessary columns from the ages dataframe to each df in the list
  LD_positions_selection <- lapply(LD_positions_sign, function(df) {
    left_join(df, ages_data[, c("SNP", "P_X", "FDR")], by = c("SNP")
    )
  })
  
  #Add the necessary columns from the gwas summary stats to each df in the list
  LD_positions_selection_gwas <- lapply(LD_positions_selection, function(df) {
    left_join(df, gwas[, c("SNP", "beta", "pval")], by = c("SNP")
    )
  })
  
  #Remove all clumps that have 0 rows after the filtering
  LD_positions_selection_gwas_clean <- Filter(
    function(df) nrow(df) > 0,
    LD_positions_selection_gwas
  )
  
  #Creating a new column, containing a combination of snp1 and snp2 (to increase speed)
  df_all_chr$key1 <- paste(df_all_chr$snp1, df_all_chr$snp2, sep=":")
  df_all_chr$key2 <- paste(df_all_chr$snp2, df_all_chr$snp1, sep=":")
  
  #Create a lookup vector, containing the r^2 values and keys (to increase speed)
  r2_lookup  <- setNames(
    c(df_all_chr$r2, df_all_chr$r2),
    c(df_all_chr$key1, df_all_chr$key2))
  
  #Use lapply to loop over the dataframes in the list and get the r^2 value between each of the SNPs in the df and the lead SNP
  LD_positions_selection_gwas_ld <- lapply(seq_along(LD_positions_selection_gwas_clean), function(i) {
    
    df <- LD_positions_selection_gwas_clean[[i]]
    
    lead_snp <- names(LD_positions_selection_gwas_clean)[i]
    
    cat("Processing", i, "of", length(LD_positions_selection_gwas_clean), ":", lead_snp, "\n")
    
    df$LD_with_lead_snp <- r2_lookup[paste(lead_snp, df$rsid, sep=":")]
    
    df
  })
  
  #Sorting dataframes in the list by pval
  LD_positions_selection_gwas_ld_pval <- lapply(LD_positions_selection_gwas_ld, function(df) {
    df[order(df$pval), ]
  })
  
  #Saving object for later
  saveRDS(
    LD_positions_selection_gwas_ld_pval,
    paste0("LD_positions_annotated_", trait, ".rds")
  )
  
}