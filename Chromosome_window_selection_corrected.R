library(data.table)
library(dplyr)
library(GenomicRanges)

urate_gwas_sign <- urate_gwas_separate_columns %>%
  filter(pval < 5e-8)

chr_lengths <- c("1"= 249250621, "2"= 243199373, "3"= 198022430, "4"= 191154276, "5"= 180915260, "6"= 171115067, "7"= 159138663, "8"= 146364022,
                 "9"= 141213431, "10"= 135534747, "11"= 135006516, "12"= 133851895, "13"= 115169878, "14"= 107349540, "15"= 102531392, "16"= 90354753,
                 "17"= 81195210, "18"= 78077248, "19"= 59128983, "20"= 63025520, "21"= 48129895, "22"= 51304566)

chr1 <- subset(urate_gwas_sign, chromosome == 1)
chr2 <- subset(urate_gwas_sign, chromosome == 2)
chr3 <- subset(urate_gwas_sign, chromosome == 3)
chr4 <- subset(urate_gwas_sign, chromosome == 4)
chr5 <- subset(urate_gwas_sign, chromosome == 5)
chr6 <- subset(urate_gwas_sign, chromosome == 6)
chr7 <- subset(urate_gwas_sign, chromosome == 7)
chr8 <- subset(urate_gwas_sign, chromosome == 8)
chr9 <- subset(urate_gwas_sign, chromosome == 9)
chr10 <- subset(urate_gwas_sign, chromosome == 10)
chr11 <- subset(urate_gwas_sign, chromosome == 11)
chr12 <- subset(urate_gwas_sign, chromosome == 12)
chr13 <- subset(urate_gwas_sign, chromosome == 13)
chr14 <- subset(urate_gwas_sign, chromosome == 14)
chr15 <- subset(urate_gwas_sign, chromosome == 15)
chr16 <- subset(urate_gwas_sign, chromosome == 16)
chr17 <- subset(urate_gwas_sign, chromosome == 17)
chr18 <- subset(urate_gwas_sign, chromosome == 18)
chr19 <- subset(urate_gwas_sign, chromosome == 19)
chr20 <- subset(urate_gwas_sign, chromosome == 20)
chr21 <- subset(urate_gwas_sign, chromosome == 21)
chr22 <- subset(urate_gwas_sign, chromosome == 22)

chromosomes <- list(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22)

merged_gwas <- list()

for (i in 1:length(chromosomes)) {
  
  chr_data <- chromosomes[[i]]
  
  windows <- unlist(slidingWindows(
    GRanges(seqnames = paste0("chr", i),
            ranges = IRanges(1, chr_lengths[i])),
    width = 3000001,
    step = 1000000
  ))
  
  mcols(windows)$window_id <- paste0("chr", i, "_", seq_along(windows))
  
  g_ranges <- GRanges(
    seqnames = paste0("chr", i),
    ranges = IRanges(start = chr_data$position, end = chr_data$position),
    SNP = chr_data$SNP,
    P = chr_data$pval
  )
  
  hits <- findOverlaps(g_ranges, windows)
  
  merged_chr <- data.frame(
    SNP = mcols(g_ranges)$SNP[queryHits(hits)],
    P = mcols(g_ranges)$P[queryHits(hits)],
    position = start(g_ranges)[queryHits(hits)],
    window_id = mcols(windows)$window_id[subjectHits(hits)],
    start = start(windows)[subjectHits(hits)],
    end = end(windows)[subjectHits(hits)],
    chromosome = rep(i, length(queryHits(hits)))
  )
  
  merged_gwas[[i]] <- merged_chr
}

merged_gwas_all <- bind_rows(merged_gwas)

merged_gwas_all$filename_npz <- paste0("http://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/", "chr", merged_gwas_all$chromosome, "_",
                                       merged_gwas_all$start, "_",
                                       merged_gwas_all$end, ".npz"
)

files_to_download_npz <- unique(merged_gwas_all$filename_npz)