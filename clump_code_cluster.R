## Retrieve the GWAS summary statistics files for the traits of interest

## Example for the trait 'urate'
wget="wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30880_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz -O 30880_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz"

## function to perform clumping of GWAS loci
clump.function <- function(wget,chrs,cores=1,ld.threshold=0.1){
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(parallel)
    library(reticulate)
    
    if(!file.exists(strsplit(wget,"-O ")[[1]][2])) system(wget)
    
    ## read in and modify GWAS file
    gwas <- fread(cmd=paste0("zcat ",strsplit(wget,"-O ")[[1]][2]),data.table=F)
    gwas <- gwas %>% separate(variant, into = c("chromosome", "position", "reference", "alternative"), sep = ":")
    gwas$SNP <- paste(gwas$chromosome, gwas$position, gwas$reference, gwas$alternative, sep = ":")
    
    ## write.csv(urate_gwas_sep, "urate_gwas_separate_columns.csv")

    gwas <- gwas %>%
    mutate(chromosome = as.numeric(chromosome)) %>%
    filter(pval < 5e-8) %>%
    filter(chromosome %in% c(1:22))

    ## all gz files
    ## gz <- mclapply(dir()[grep("^chr.*\\.gz$",dir())],function(i) fread(i,data.table=F),mc.cores=cores)
    ## save(gz,file="gz_file_info.rdat")

    load("gz_file_info.rdat")

    ## find gz files overlapping significant GWAS SNPs
    gz.sig <- dir()[grep("^chr.*\\.gz$",dir())][sapply(gz,function(x) sum(paste0(x[,2],":",x[,3],":",x[,4],":",x[,5]) %in% paste0(gwas[,1],":",gwas[,2],":",gwas[,3],":",gwas[,4])),USE.NAMES=F)>0]
    
    np <- import("numpy")

    ##Create the necessary empty lists
    df_all_chr <- vector("list")
    meta_all_chr <- vector("list")
        
    for (chr in chrs) {
        cat("analyzing chromosome ",chr,"\n")
                
        ##Create a vector containing the names of the gz and npz files respectively
        gz_files <- grep(paste0("^chr", chr, "_.*\\.gz$"), gz.sig, value = TRUE)

        download <- mclapply(gz_files, function(i) {
            base <- sub("\\.gz$", "", i)
            
            cmd <- paste0(
                "aws s3 cp s3://broad-alkesgroup-ukbb-ld/UKBB_LD/ . ",
                "--recursive ",
                "--exclude \"*\" ",
                "--include \"", base, ".npz\" ",
                "--include \"", base, ".npz2\" ",
                "--no-sign-request"
            )
            
            system(cmd)
        },mc.cores=20)
           
        npz_files <- grep(paste0("^chr", chr, "_.*\\.npz(2)?$"), dir("."), value = TRUE)
        
        final_df <- vector("list")
        meta_df <- vector("list")

        for (i in 1:length(gz_files)) {
            cat("npz file:",i," (of ",length(gz_files),")\r")
            npz_file <- npz_files[i]
            gz_file <- gz_files[i]

            ldfile <- np$load(npz_file)
            rows <- as.integer(ldfile$f[["row"]])
            cols <- as.integer(ldfile$f[["col"]])
            vals <- as.numeric(ldfile$f[["data"]])
            vals_sq <- vals^2

            meta <- fread(gz_file,
                          sep = "\t",
                          data.table = FALSE
                          )

            snp_ids <- meta$rsid

            df <- data.frame(
                snp1 = snp_ids[rows + 1],
                snp2 = snp_ids[cols + 1],
                pos1 = meta$position[rows + 1],
                pos2 = meta$position[cols + 1],
                r2   = vals_sq
            )

            ld_threshold <- ld.threshold

            df_high_LD <- subset(df, r2 > ld_threshold)
            df_sign <- df_high_LD[df_high_LD$pos1 %in% gwas$position[gwas$chromosome==chr] & df_high_LD$pos2 %in% gwas$position[gwas$chromosome==chr], ]

            final_df[[i]] <- df_sign
            meta_df[[i]] <- meta
        }
        
        df_total <- bind_rows(final_df)
        df_total_unique <- unique(df_total)
        df_total_unique <- df_total_unique[df_total_unique[,1]!=df_total_unique[,2],]
        
        meta_total <- bind_rows(meta_df)
        meta_total_unique <- unique(meta_total)        
        
        df_all_chr[[chr]] <- df_total_unique
        meta_all_chr[[chr]] <- meta_total_unique

        save(df_total_unique,meta_total_unique,file=paste0(strsplit(wget,"-O ")[[1]][2],"_chr",chr,"_df.rdat"))
        
        ## remove npz files for a given chromosome
        remove <- sapply(npz_files,function(x) system(paste0("rm ",x)))
    }
    
    df_all_chr <- bind_rows(df_all_chr)
    meta_all_chr <- bind_rows(meta_all_chr)    
    
    ## save LD information
    save(df_all_chr,file=paste0(strsplit(wget,"-O ")[[1]][2],".rdat"))
    
    ## clumping algorithm
    ## order GWAS file by p value
    gwas <- gwas[order(gwas$pval), ]
    
    gwas_annot <- merge(
        gwas,
        meta_all_chr,
        by = c("chromosome", "position")
    )
    
    gwas_annot_sorted <- gwas_annot[order(gwas_annot$pval), ]
    
    gwas_annot_sorted_b <- gwas_annot[order(abs(gwas_annot$beta), decreasing = TRUE), ]
    
    LD_sets <- vector("list")
    LD_positions <- vector("list")
    
    step=1
    
    while (nrow(gwas_annot_sorted) > 0) {
        cat("clump number: ",step,"\n")
        step = step + 1
        
        lead_snp <- gwas_annot_sorted$rsid[1]
        lead_pval <- gwas_annot_sorted$pval[1]
        
        if (lead_pval == 0) {
            lead_snp <- gwas_annot_sorted_b$rsid[gwas_annot_sorted_b$pval==0][1]
            lead_pval <- gwas_annot_sorted_b$pval[gwas_annot_sorted_b$pval==0][1]
        }
        
        snps_in_set <- lead_snp
        previous_size <- 0
        sig_snps <- gwas_annot$rsid
        
        step2 = 1
        
        while (length(snps_in_set) > previous_size) {
            cat("LD clump iteration: ",step2,"\n")
            step2 = step2 + 1
            
            previous_size <- length(snps_in_set)
            
            new_snps <- unique(unlist(c(
                df_all_chr[df_all_chr$snp1 %in% snps_in_set | df_all_chr$snp2 %in% snps_in_set,c("snp1","snp2")]
            )))
            
            snps_in_set <- unique(c(snps_in_set, new_snps))
            
        }
        
        LD_sets[[lead_snp]] <- snps_in_set
        LD_positions[[lead_snp]] <- meta_all_chr[match(snps_in_set, meta_all_chr$rsid),]
        
        gwas_annot_sorted <- gwas_annot_sorted[!(gwas_annot_sorted$rsid %in% snps_in_set), ]
        gwas_annot_sorted_b <- gwas_annot_sorted_b[!(gwas_annot_sorted_b$rsid %in% snps_in_set), ]

    }
    ## save data objects
    save(LD_sets,LD_positions,file=paste0(strsplit(wget,"-O ")[[1]][2],"_clumps.rdat"))
    
    ## remove GWAS file
    system(paste0("rm ",strsplit(wget,"-O ")[[1]][2]))
}
