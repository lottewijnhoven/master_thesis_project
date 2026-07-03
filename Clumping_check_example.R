#Load the necessary packages
library(dplyr)
library(tidyr)
library(data.table)
library(dplyr)
library(qqman)

#Make sure there is a list available, called 'LD_positions', that contains the clumps including their information
#namely: rsid, chromosome, position, allele 1 and allele 2, where the name of each clump is the respective lead SNP

#Read in the GWAS summary statistics for the trait of interest (in this case 'urate')
gwas <- fread("30880_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz")

#Separate the 'variant' column into 'chromosome', 'position', 'reference' and 'alternative' columns
gwas <- gwas %>% separate(variant, into = c("chromosome", "position", "reference", "alternative"), sep = ":")

#Recreate the original 'variant' column, and call it 'SNP'
gwas$SNP <- paste(gwas$chromosome, gwas$position, gwas$reference, gwas$alternative, sep = ":")

#Create a subset of the summary statistics for the chromosome of interest (in this case '22')
gwas_chr <- subset(gwas, chromosome == 22)

#Indicate the length (in bp) of the chromosome of interest
chr_length <- 115169878

#Ensure that the position is considered a numeric variable
gwas_chr$position <- as.numeric(gwas_chr$position)
gwas_chr$chromosome <- as.numeric(gwas_chr$chromosome)

#Give the pvalues noted down as '0' an extremely low value instead
gwas_chr$pval[gwas_chr$pval == 0] <- .Machine$double.xmin

#Create a new 'x_coordinate' column, containing the chromosomal positions
gwas_chr$x_coordinate <- gwas_chr$position

#Check if any NA's or infinite values are still present
summary(gwas_chr$pval)
anyNA(gwas_chr$pval)          
any(!is.finite(gwas_chr$pval)) 

#Rename the chromosome, position and pval columns and create a new object containing only these three columns
gwas_qqman <- gwas_chr %>%
  rename(CHR = chromosome, BP = position, P = pval) %>%
  select(CHR, BP, P)

gwas_qqman$CHR <- as.numeric(gwas_qqman$CHR)

#Add a SNP column to the df
gwas_qqman$SNP <- gwas_chr$SNP

#Plot the manhattan plot (in this case for urate, chromosome 22)
manhattan(gwas_qqman,
          chr = "CHR",
          bp = "BP", 
          p = "P",
          snp = "SNP",
          main = "Manhattan Plot of urate GWAS for chr22",
          col = c("darkgrey", "black"),
          cex = 1,
          genomewideline = -log10(5e-8),
          suggestiveline = FALSE,
          xlim = c(min(gwas_qqman$BP),
                   max(gwas_qqman$BP)),
          ylim = c(0, 18))

#Keep only the clumps in the 'LD_positions' list that are located on the chromosome of interest (in this case 22)

keep <- sapply(LD_positions, function(df) 22 %in% df$chromosome)
chr22_clumps <- LD_positions[keep]

#Create a copy of the generated list, but this time only containin the positions of the SNPs
positions_list <- lapply(chr22_clumps, function(df) as.character(df$position))

#Determine the colours used to highlight the clumps in the Manhattan plot later
cols <- rainbow(length(positions_list))

#Loop over the different clumps for the chromosome of interest, get the necessary data from the GWAS summary statistics
#and plot the SNPs belonging to each of the clumps on top of the generated Manhattan plot (in one colour per clump)
for (i in seq_along(positions_list)) {
  clump_data <- gwas_qqman %>%
    filter(BP %in% positions_list[[i]])
  
  points(
    clump_data$BP,
    -log10(clump_data$P),
    col = cols[i],
    pch = 19,
    cex = 1.5
  )
}


