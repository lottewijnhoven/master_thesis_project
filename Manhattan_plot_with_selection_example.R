#Load the necessary libraries
library(dplyr)
library(tidyr)
library(data.table)
library(dplyr)
library(ggplot2)

#Read in the GWAS summary statistics for the trait of interest (in this case 'urate')
gwas <- fread("30880_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz")

#Load the AGES summary statistics for later
ages_data <- fread("Selection_Summary_Statistics_01OCT2025.tsv.gz")
ages_data$SNP <- paste(ages_data$CHROM, ages_data$POS, ages_data$REF, ages_data$ALT, sep = ":")

#Separate the 'variant' column into 'chromosome', 'position', 'reference' and 'alternative' columns
gwas <- gwas %>% separate(variant, into = c("chromosome", "position", "reference", "alternative"), sep = ":")

#Recreate the original 'variant' column, and call it 'SNP'
gwas$SNP <- paste(gwas$chromosome, gwas$position, gwas$reference, gwas$alternative, sep = ":")

#Create a vector containing the lengths of all autosomes
chr_lengths <- c("1"= 249250621, "2"= 243199373, "3"= 198022430, "4"= 191154276, "5"= 180915260, "6"= 171115067, "7"= 159138663, "8"= 146364022,
                 "9"= 141213431, "10"= 135534747, "11"= 135006516, "12"= 133851895, "13"= 115169878, "14"= 107349540, "15"= 102531392, "16"= 90354753,
                 "17"= 81195210, "18"= 78077248, "19"= 59128983, "20"= 63025520, "21"= 48129895, "22"= 51304566)

#Create a vector containing the offsets of each of the chromosomes on the x axis
chr_offsets <- cumsum(chr_lengths) - chr_lengths

#Make sure the position column is numerical
gwas$position <- as.numeric(gwas$position)

#Give the pvalues noted down as '0' an extremely low value instead
gwas$pval[gwas$pval == 0] <- .Machine$double.xmin

#Filter to include only GWAS summary statistics for the autosomes
gwas <- gwas %>%
  filter(chromosome %in% names(chr_lengths))

#Create a new 'x_coordinate' column, containing the x coordinates (chromosomal positions) of all SNPs
gwas$x_coordinate <- 
  gwas$position + 
  chr_offsets[as.character(gwas$chromosome)]

#Add the necessary columns from the AGES summary statistics (P_X and FDR) to the GWAS data
gwas_ages <- left_join(gwas, ages_data[, c("SNP", "P_X", "FDR")], by = c("SNP"))

#Add a column to the data where SNPs are grouped into three categories, based on selection (for colouring later on)
gwas_ggplot <- gwas_ages %>%
  mutate(group = case_when(
    P_X < 0.05 & FDR < 0.05 ~ "Significant",
    P_X < 0.05 & FDR >= 0.05 ~ "P<0.05, FDR>0.05",
    P_X >= 0.05 & FDR >= 0.05 ~ "Not significant",
    TRUE ~ "Not in AGES"
  ))

#Determine the locations of the centres of the chromosomes on the x axis
chr_centers <- chr_offsets + chr_lengths / 2

#Create vectors for each of the categories, containing the SNPs that belong to the category (example here for urate)
GxE <- "4:103112470:A:G"
vQTL <- "2:27730940:T:C"
both <- c("4:89054667:A:G", "4:89052323:G:T")
no_sel_GxE <- "4:9936437:G:C"
no_sel_vQTL <- "6:25759066:G:T"
not_in_AGES_GxE <- "6:81433192:G:T"
not_in_AGES_vQTL <- "4:9970570:A:G"

#For each of SNPs in the categories, get the necessary data from the dataset
snps_GxE <- gwas_ggplot %>% filter(SNP %in% GxE)
snps_vQTL <- gwas_ggplot %>% filter(SNP %in% vQTL)
snps_both <- gwas_ggplot %>% filter(SNP %in% both)
snps_ns_GxE <- gwas_ggplot %>% filter(SNP %in% no_sel_GxE)
snps_ns_vQTL <- gwas_ggplot %>% filter(SNP %in% no_sel_vQTL)
snps_no_AGES_GxE <- gwas_ggplot %>% filter(SNP %in% not_in_AGES_GxE)
snps_no_AGES_vQTL <- gwas_ggplot %>% filter(SNP %in% not_in_AGES_vQTL)

#Give each of the SNPs in the categories a corresponding label (in a new column, named 'SNP_set')
gwas_ggplot <- gwas_ggplot %>%
  mutate(SNP_set = case_when(
    SNP %in% GxE ~ "GxE",
    SNP %in% vQTL ~ "vQTL",
    SNP %in% both ~ "Both",
    SNP %in% no_sel_GxE ~ "No sel GxE",
    SNP %in% no_sel_vQTL ~ "No sel vQTL",
    SNP %in% not_in_AGES_GxE ~ "Not in AGES GxE",
    SNP %in% not_in_AGES_vQTL ~ "Not in AGES vQTL",
    TRUE ~ "Other"
  ))

#Give each SNP a number, corresponding to the order in which the datapoints will be plotted on top of the graph
gwas_ggplot <- gwas_ggplot %>%
  mutate(plot_order = case_when(
    group == "Significant" ~ 3,
    group == "P<0.05, FDR>0.05" ~ 2,
    group == "Not significant" ~ 1,
    TRUE ~ 0
  )) %>%
  arrange(plot_order) #And sort the data by this variable

#Create the Manhattan plot, and plot the SNPs of interest on top in the colors and shapes corresponding to their categories
ggplot(gwas_ggplot, aes(x = x_coordinate, y = -log10(pval), color = group)) +
  geom_point(size = 1) +
  geom_point(
    data = subset(gwas_ggplot, SNP_set != "Other"),
    aes(shape = SNP_set, color = group),
    size = 4, stroke = 1,
    show.legend = TRUE
  ) +
  scale_x_continuous(
    breaks = chr_centers,
    labels = names(chr_lengths)
  ) +
  scale_color_manual(values = c(
    "Significant" = "darkgreen",
    "P<0.05, FDR>0.05" = "lightblue",
    "Not significant" = "pink",
    "Not in AGES" = "grey"
  )) +
  scale_shape_manual(values = c(
    "GxE" = 17, "vQTL" = 15, "Both" = 18,
    "No sel GxE" = 3, "No sel vQTL" = 7,
    "Not in AGES GxE" = 8, "Not in AGES vQTL" = 4
  )) +
  labs(
    title = "Manhattan plot of urate GWAS with selection",
    x = "Chromosome",
    y = "-log10(P-value)",
    shape = "Highlighted SNPs",
    color = "Significance"
  ) +
  theme_minimal()
  
