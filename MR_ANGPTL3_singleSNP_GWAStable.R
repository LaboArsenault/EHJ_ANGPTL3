#!/usr/bin/env Rscript

################################################################################
#i. Script options ####
library(optparse)
{
  option_list = list(
## GLOBAL OPTIONS
    make_option("--wd",action="store",default="./",type="character",
                help="working directory from which to access files"),
    make_option("--window", action="store", default=5e05, type='numeric',
                help=" inclusion and clumping window (default 500 kb)"),
    make_option("--pval", action="store", default=5e-08, type='numeric',
                help=" pval threshold (default 5e-08)"),
    make_option("--r2", action="store", default=0.1, type='numeric',
                help=" clumping threshold (default 0.1)"),
    make_option("--from", action="store", default=NA, type='numeric',
                help=" If you want to run this script in multiple sequences, select which min protein number to include"),
    make_option("--to", action="store", default=NA, type='numeric',
                help=" If you want to run this script in multiple sequences, select which max protein number to include"),
    make_option("--study", action="store", default="FinnGen_R7", type='character',
                help=" which study/consortium as outcome ? (UKB (Zhou) or FinnGen_R7 (default))"),
    make_option("--gwas_table", action="store", default="/mnt/sda/gobemi01/Data_instruments/gwas_table_07-10-2022.xlsx", type='character',
                help=" study directory where GWAS_table file is stored"),
    make_option("--genes_ref_file", action="store", default="./data/GRCh38_chr_pos_genes.txt", type='character',
                help=" path to genes reference file (name, pos, chr) ( HGNC Build 38 (default))"),
    make_option("--plink_binaries", action="store", default="./data/EUR_rs", type='character',
                help=" (Optional) path to plink binary files for local clumping "),
    make_option("--r2_corr", action="store", default=0.6, type='numeric',
                help=" (Optional) R2 threshold for analyses including LD-correction (default 0.6) "),
    make_option("--proxies", action="store", default=TRUE, type='logical',
                help=" (Optional) Should you use proxies ? (default TRUE)"),
    make_option("--reverse", action="store", default=FALSE, type='logical',
                help=" (Optional) Should you perform reverse PheWAS MR ? (default FALSE)"),
    make_option("--cis_trans", action="store", default="cis", type='character',
                help=" Do you want only cis or cis_trans SNP as exposure? cis (default cis)"),
    make_option("--singleSNP_multiSNP", action="store", default="singleSNP", type='character',
                help=" Perform single SNP or multi SNP phenome-wide MR analysis? (default singleSNP)"),

    
## EXPOSURE OPTIONS
  make_option("--split_var", action="store", default="exposure", type='character',
              help=" The colname to split the all_out file by (default exposure)"),
  make_option("--exp_gwas_file", action="store", default=NA, type='character',
            help="File for exposure(s)"),
  make_option("--exp_pheno", action="store", default="ANGPTL3", type='character',
            help=" Exposure phenotype name (default ANGPTL3)"),
make_option("--exp_start_pos", action="store", default=1, type='integer',
            help="Gene start pos_b38"),
make_option("--exp_end_pos", action="store", default=1000000, type='integer',
            help="Gene end pos_b38"),
  make_option("--exp_beta_col", action="store", default="Beta", type='character',
            help=" beta col"),
  make_option("--exp_se_col", action="store", default="SE", type='character',
            help=" SE col"),
  make_option("--exp_chr_col", action="store", default=NA, type='character',
            help=" chromosome col"),
  make_option("--exp_pos_col", action="store", default=NA, type='character',
            help=" position col"),
  make_option("--exp_pval_col", action="store", default="Pval", type='character',
            help=" P-value col"),
  make_option("--exp_ea_col", action="store", default="effectAllele", type='character',
            help=" effect allele col"),
  make_option("--exp_oa_col", action="store", default="otherAllele", type='character',
            help=" other allele col"),
  make_option("--exp_eaf_col", action="store", default=NA, type='character',
            help=" effect allele frequency col"),
  make_option("--exp_ss_col", action="store", default=NA, type='character',
            help=" sample size col"),
  make_option("--exp_snp_col", action="store", default="SNP", type='character',
            help=" SNP rsids column name (default SNP)"),
  make_option("--exp_ncase", action="store", default=NA, type='numeric',
            help=" Cases number for exposure"),
  make_option("--exp_ncontrol", action="store", default=NA, type='numeric',
            help=" Controls number for exposure"),
  make_option("--exp_ss_total", action="store", default=NA, type='numeric',
            help=" Total samplesize for exposure"),
  make_option("--exp_p_is_log", action="store", default=FALSE, type='logical',
            help=" (Optional) Is exposure pvalue in log10 ? "),
  make_option("--exp_is_vcf", action="store", default=FALSE, type='logical',
            help=" is exposure file a vcf file ? "),
  make_option("--exp_vcf_skip", action="store", default=1, type='numeric',
            help=" number of header lines to skip in vcf file. Do not specify anything if file is not vcf, otherwise CoJo will crash. "),
  make_option("--exp_snp_col_num", action="store", default=3, type='numeric',
            help=" exposure SNP column number "),
  make_option("--exp_chr_col_num", action="store", default=1, type='numeric',
            help=" exposure chromosome column number "),
  make_option("--exp_pos_col_num", action="store", default=2, type='numeric',
            help=" exposure position column number "),
  make_option("--exp_a1_col_num", action="store", default=4, type='numeric',
            help=" exposure allele 1 column number "),
  make_option("--exp_a2_col_num", action="store", default=5, type='numeric',
            help=" exposure allele 2 column number "),
  make_option("--exp_b_col_num", action="store", default=12, type='numeric',
            help=" exposure beta column number "),
## IF FILE IS VCF, the next five columns number should be specified
## according to their position in the FORMAT column, not in
  make_option("--exp_maf_col_num", action="store", default=4, type='numeric',
            help=" exposure maf column number "),
  make_option("--exp_beta_col_num", action="store", default=5, type='numeric',
            help=" exposure beta column number "),
  make_option("--exp_se_col_num", action="store", default=6, type='numeric',
            help=" exposure se column number "),
  make_option("--exp_p_col_num", action="store", default=7, type='numeric',
            help=" exposure pval column number "),
  make_option("--exp_n_col_num", action="store", default=8, type='numeric',
            help=" exposure N column number "),


## OUTCOME OPTIONS
  make_option("--min_ratio", action="store", default="0.001", type='numeric',
            help=" Minimum ratio case/control to keep? (default 0.001)"),
  make_option("--min_cases", action="store", default="500", type='numeric',
            help=" Minimum ratio case/control to keep? (default 500)"),
  make_option("--to_exclude_file", action="store", default=NA, type='character',
            help=" Path to to_exclude.txt file if we want to remove phenotype from analysis? (default NA)"),
  make_option("--out_is_vcf", action="store", default=FALSE, type='logical',
            help=" is exposure file a vcf file ? "),
  make_option("--out_vcf_skip", action="store", default=NA, type='numeric',
            help=" number of header lines to skip in vcf file. Do not specify anything if file is not vcf. "),
  make_option("--out_beta_col", action="store", default="Beta", type='character',
            help=" beta col"),
  make_option("--out_se_col", action="store", default=NA, type='character',
            help=" SE col"),
  make_option("--out_chr_col", action="store", default=NA, type='character',
            help=" chromosome col"),
  make_option("--out_pos_col", action="store", default=NA, type='character',
            help=" position col"),
  make_option("--out_pval_col", action="store", default=NA, type='character',
            help=" P-value col"),
  make_option("--out_ea_col", action="store", default=NA, type='character',
            help=" effect allele col"),
  make_option("--out_oa_col", action="store", default=NA, type='character',
            help=" other allele col"),
  make_option("--out_eaf_col", action="store", default=NA, type='character',
            help=" effect allele frequency col"),
  make_option("--out_ss_col", action="store", default=NA, type='character',
            help=" sample size col"),
  make_option("--out_snp_col", action="store", default="rsids", type='character',
            help=" SNP rsids col"),
  make_option("--out_p_is_log", action="store", default=FALSE, type='logical',
            help=" (Optional) Is outcome pvalue in log10 ? "),
  make_option("--out_ss_total", action="store", default=NA, type='numeric',
            help=" Total samplesize for outcome"),


## OUTPUT OPTIONS
  make_option("--save_every", action="store", default=1, type='numeric',
            help=" save results every X analysis (default 1)"),
  make_option("--outname_res", action="store", default=NA, type='character',
            help=" output name of results"),
  make_option("--outname_dat", action="store", default=NA, type='character',
            help=" output name of MR data (harmonized / clumped)"),
  make_option("--outname_coloc", action="store", default=NA, type='character',
             help=" output name of coloc results"),
  make_option("--save_dat", action="store", default=TRUE, type='logical',
            help=" should you save MR data for each protein? (default TRUE) "),
  make_option("--check_memory", action="store", default=TRUE, type='logical',
            help=" (Optional) Should you check for memory usage and kill the task if it exceeds a certain amount of memory ? (default TRUE)"),
  make_option("--check_memory_lim", action="store", default=0.9, type='numeric',
            help=" (Optional)  If --check_memory, what is the ratio limit (current_mem_use/total_mem) before killing the task (Maximum is 0.95) ? (default 0.9)"),
  make_option("--nthreads", action="store", default=1, type='numeric',
            help=" How many threads should you use to perform analyses ? (default 1) ")   
  )
  opt = parse_args(OptionParser(option_list=option_list))   
}
############################################################################## #
############################################################################## #
#### ii. Loading libraries ####
library(data.table)
library(TwoSampleMR)
library(MRInstruments)
library(xlsx)
library(stringr)
library(rJava)
library(PheWAS)
library(ggplot2)
library(ieugwasr)
library(coloc)
library(MendelianRandomization)
setwd(opt$wd)
setDTthreads(threads = opt$nthreads, throttle = 1024*opt$nthreads)
############################################################################## #
############################################################################## #
#### 1. Parameters ####
{
  for(i in 1:length(opt)){
    if(!is.na(opt[[i]])){
      if(opt[[i]] == "NA"){
        opt[[i]] <- NA
      }
    }
  }
  from_na <- !is.na(opt$from)
  to_na <- !is.na(opt$to)
  use_from_to <- ifelse(test = (from_na+to_na)!=1,
                        yes = TRUE,
                        no = FALSE)
  if(!use_from_to){
    stop(" You specified only one of the 'from' and 'to' arguments. Specify either none or both of them.")
  }
  if(!(grepl(pattern = "UKB", x = opt$study) | grepl(pattern = "FinnGen", x = opt$study) )){
    message(" Study name provided is not UKB or FinnGen. This script has not been tested yet on other studies and may be unstable for this analysis. ")
  }
  if(is.na(opt$exp_ss_total)){
    stop("Total samplesize for exposure is required.")
  }
  window <- as.numeric(opt$window)
  pval_threshold <- as.numeric(opt$pval)
  clump_r2 <- as.numeric(opt$r2)
  suff <- paste0(opt$cis_trans,"_",opt$singleSNP_multiSNP)
  if(is.na(opt$outname_res)){
    if(use_from_to & (from_na+to_na==2)){
      outname <- paste0("./results/MRres/", opt$study,"/",opt$exp_pheno, "/w",opt$window/1000,"_r2_",opt$r2,"/MR_results_",opt$exp_pheno,"_",suff,"_pval.", pval_threshold,"_w",opt$window/1000, "_r2.", clump_r2, "_from.", opt$from, "_to.", opt$to, ".txt")
    } else {
      outname <- paste0("./results/MRres/", opt$study,"/",opt$exp_pheno,"/w",opt$window/1000,"_r2_",opt$r2,"/MR_results_",opt$exp_pheno,"_",suff,"_pval.", pval_threshold,"_w",opt$window/1000,"_r2.", clump_r2, ".txt")
    }
  } else {
    outname <- opt$outname_res
  }
  if(opt$save_dat){
    if(is.na(opt$outname_dat)){
      outname_dat <- paste0("./results/MRdat/", opt$study, "/",opt$exp_pheno,"/w",opt$window/1000,"_r2_",opt$r2,"/MR_dat_")
    } else {
      outname_dat <- opt$outname_dat
    }
  }
  if(is.na(opt$outname_coloc)){
    outname_coloc <- paste0("./results/coloc/", opt$study,"/",opt$exp_pheno,"/w",opt$window/1000,"_r2_",opt$r2, "/coloc_res_")
  } else {
    outname_coloc <- opt$outname_coloc
  }
  
  for(directory in c("MRres", "MRdat", "coloc")){
    if(!dir.exists(paste0("./results/", directory, "/", opt$study,"/",opt$exp_pheno, "/" ,"w",opt$window/1000,"_r2_",opt$r2,"/"))){
      dir.create(path = paste0("./results/", directory, "/", opt$study,"/",opt$exp_pheno, "/","w",opt$window/1000,"_r2_",opt$r2, "/"))
    }
  }
}
############################################################################## #
############################################################################## #
#### 2. GWAS_table directory ####
{
  gwas_table <- readxl::read_xlsx(opt$gwas_table)
}
############################################################################## #
#### 3. Exclusion and Selection of Phenotypes ####
{
  min_ratio <- as.numeric(opt$min_ratio)
  #gwas_table <- gwas_table[which(gwas_table$ratio >= min_ratio),]
  if(is.na(opt$to_exclude_file)){
    phewas_dataset <- gwas_table
  } else {
    to_exclude <- as.data.frame(fread(opt$to_exclude_file, header = TRUE, stringsAsFactors = FALSE,nThread = opt$nthreads))
    phewas_dataset <- gwas_table[(-which(gwas_table$outcome %in% to_exclude$outcome)),]
  }
}
############################################################################## #
############################################################################## #
#### 4. Instrument (exposure) GWAS ####
{
  if(grepl(pattern = "deCODE", x = opt$exp_pheno)){
   exp_gwas <- as.data.frame(fread(file = opt$exp_gwas_file, header = TRUE, stringsAsFactors = FALSE, nThread = opt$nthreads))
  exp_gwas$phenotype <- opt$exp_pheno
  exp_gwas$samplesize <- opt$exp_ss_total
  exp_gwas$Chrom <- str_extract_all(string=exp_gwas[,opt$exp_chr_col],pattern="\\d+",simplify = TRUE)
  exp_gwas$Chrom <- as.numeric(exp_gwas$Chrom)
  exp_gwas <- dplyr::arrange(exp_gwas,exp_gwas[,opt$exp_pval_col])
  exp_gwas$Pval <- as.numeric(exp_gwas[,opt$exp_pval_col])
  exp_dat <- exp_gwas[which(exp_gwas$Pval <= as.numeric(opt$pval)),]
  # Format GWAS for COLOC
  exp_gwas <- exp_gwas[which(exp_gwas[,opt$exp_chr_col] == exp_gwas[1,opt$exp_chr_col]),]
  exp_gwas <- exp_gwas[which((exp_gwas[,opt$exp_pos_col] >= (as.integer(opt$exp_start_pos)-opt$window/2)) & (exp_gwas[,opt$exp_pos_col] <= as.integer(opt$exp_end_pos)+opt$window/2)),]
  exp_gwas <- subset(x = exp_gwas, select = c("rsids", "Beta", "SE", "effectAllele", "otherAllele", "Chrom", "Pos", "Pval","ImpMAF"))
  colnames(exp_gwas) <- c("SNP", "beta", "se", "ea", "oa", "chr", "pos", "pval","maf")
  exp_gwas <- na.omit(exp_gwas) %>% unique(.) %>% .[order(.$pval),] %>% filter(!duplicated(SNP))
  exp_gwas$pheno <- paste0(opt$exp_pheno)
  exp_gwas$SS <- 35559
  exp_gwas <- format_data(exp_gwas, 
                        type = "exposure",
                        phenotype_col = "pheno",
                        snp_col = "SNP",
                        beta_col = "beta",
                        se_col = "se",
                        pval_col = "pval",
                        eaf_col = "maf",
                        effect_allele_col = "ea",
                        other_allele_col = "oa",
                        samplesize_col = "SS",
                        chr_col = "chr",
                        pos_col = "pos",
                        min_pval=1e-300
  ) %>% as.data.table(.)
  } else {
    exp_gwas <- as.data.frame(fread(file = opt$exp_gwas_file, header = TRUE, stringsAsFactors = FALSE, nThread = opt$nthreads))
    exp_gwas$phenotype <- opt$exp_pheno
    exp_gwas$samplesize <- opt$exp_ss_total
    exp_gwas$Chrom <- as.numeric(exp_gwas$chr)
    exp_gwas <- dplyr::arrange(exp_gwas,exp_gwas[,opt$exp_pval_col])
    exp_gwas$Pvalue <- as.numeric(exp_gwas[,opt$exp_pval_col])
    exp_dat <- exp_gwas[which(exp_gwas$Pvalue <= as.numeric(opt$pval)),]
    # Format GWAS for COLOC
    exp_gwas <- exp_gwas[which(exp_gwas[,opt$exp_chr_col] == exp_gwas[1,opt$exp_chr_col]),]
    exp_gwas <- exp_gwas[which((exp_gwas[,opt$exp_pos_col] >= (as.integer(opt$exp_start_pos)-opt$window/2)) & (exp_gwas[,opt$exp_pos_col] <= as.integer(opt$exp_end_pos)+opt$window/2)),]
    exp_gwas <- subset(x = exp_gwas, select = c("rsid", "Effect", "StdErr", "Allele1", "Allele2", "chr", "pos_b38", "Pvalue","Freq1"))
    colnames(exp_gwas) <- c("SNP", "beta", "se", "ea", "oa", "chr", "pos", "pval","maf")
    exp_gwas <- na.omit(exp_gwas) %>% unique(.) %>% .[order(.$pval),] %>% filter(!duplicated(SNP))
    exp_gwas$pheno <- paste0(opt$exp_pheno)
    exp_gwas$SS <- as.numeric(opt$exp_ss_total)
    exp_gwas <- format_data(exp_gwas, 
                            type = "exposure",
                            phenotype_col = "pheno",
                            snp_col = "SNP",
                            beta_col = "beta",
                            se_col = "se",
                            pval_col = "pval",
                            eaf_col = "maf",
                            effect_allele_col = "ea",
                            other_allele_col = "oa",
                            samplesize_col = "SS",
                            chr_col = "chr",
                            pos_col = "pos",
                            min_pval=1e-300
    ) %>% as.data.table(.) 
  }
  #########################
if(opt$cis_trans == "cis_selection_TG"){
  # top_snps<-c("rs1168017","rs1570694","rs11207998","rs12749263","rs34693359","rs76272805","rs74716393","rs17388017",
  #             "rs148662400","rs61775882","rs61775912","rs79152165","rs79151558","rs112548469","rs116423924","rs138783696",
  #             "rs564992513","rs115436978","rs61775945","rs12737254","rs78942743") #UKB ANGPTL3 SNPs from Xiao et al., EJPC, 2022.
  tg<-fread("/mnt/sda/gobemi01/GWAS/GWAS_in_b38/b38_logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz")
  top_snps<- tg %>% filter(pvalue <= as.numeric(5e-8))
  top_snps <- top_snps[which(as.numeric(top_snps$CHROM) == exp_dat[1,opt$exp_chr_col]),]
  top_snps <- top_snps[which((top_snps$pos_b38 >= (as.integer(opt$exp_start_pos)-opt$window/2)) & (top_snps$pos_b38 <= (as.integer(opt$exp_start_pos)+opt$window/2))),]
  exp_dat <- exp_dat %>% filter(rsids %in% top_snps$rsID)
  }
  if(opt$cis_trans == "cis"){
    exp_dat <- exp_dat[which(exp_dat[,opt$exp_chr_col] == exp_dat[1,opt$exp_chr_col]),]
    exp_dat <- exp_dat[which((exp_dat[,opt$exp_pos_col] >= (as.integer(opt$exp_start_pos)-opt$window/2)) & (exp_dat[,opt$exp_pos_col] <= as.integer(opt$exp_end_pos)+opt$window/2)),]
  }   
#### 4.1. Format exposure dat ####
  command_format <- paste0("exposure_format <- TwoSampleMR::format_data(dat = exp_dat, type = 'exposure',
                           snp_col = '", opt$exp_snp_col,
                           "', beta_col = '", opt$exp_beta_col,
                           "', se_col = '", opt$exp_se_col,
                           "', effect_allele_col = '", opt$exp_ea_col,
                           "', other_allele_col = '", opt$exp_oa_col,
                           "', pval_col = '", opt$exp_pval_col,
                           "', phenotype_col = 'phenotype' , min_pval = 1e-300" )
  if(!is.na(opt$exp_pos_col)){
    command_format <- paste0(command_format, ", pos_col = '", opt$exp_pos_col, "'")
  }
  if(!is.na(opt$exp_chr_col)){
    command_format <- paste0(command_format, ", chr_col = 'Chrom'")
  }
  if(!is.na(opt$exp_eaf_col)){
    command_format <- paste0(command_format, ", eaf_col = '", opt$exp_eaf_col, "'")
  }
  if(!is.na(opt$exp_ss_col)){
    command_format <- paste0(command_format, ", samplesize_col = '", opt$exp_ss_col, "'")
  }
  command_format <- paste0(command_format, ")")
  eval(parse(text = command_format))
  
  if(is.na(opt$exp_ss_col)){
    exposure_format$samplesize.exposure <- opt$exp_ss_total
  }
  if(opt$exp_p_is_log){
    exposure_format$pval.exposure <- 10**(exposure_format$pval.exposure)
  }
  # If performing reverse MR, keeping only SNPs with p < threshold in outcome.
  if(opt$reverse){
    exposure_all.clean <- exposure_format
    exposure_format.clean <- subset(x = exposure_format, subset = pval.exposure <= pval_threshold)
    fwrite(x = as.list(exposure_format.clean$SNP),
           file = paste0(outname, ".reverse.snps"),
           sep = "\n", append = FALSE, col.names=FALSE, row.names=FALSE)
  }
### Keeping only SNPs with maf >= 1%
exposure_format <- subset(x = exposure_format, subset = eaf.exposure >= 0.01)
############################################################################## #
############################################################################## #
#### 4.2. Clumping ####
message("\n\t Performing clumping...\n")
# Étant donné les très petites valeurs de pval qui atteignent un plancher, on se servira du abs_zscore/2 pour générer des
# pval alternatives (et plus grandes) qui serviront au clumping...
#if(opt$singleSNP_multiSNP == "multiSNP"){
  for_clumping = exposure_format
  for_clumping$good_pval=for_clumping$"pval.exposure"
  for_clumping$abs_zscore = abs(for_clumping$"beta.exposure"/for_clumping$"se.exposure")
  for_clumping$"pval.exposure"=2*pnorm(-(for_clumping$abs_zscore)/10)
  out_clump = try(ld_clump(data.frame(rsid=for_clumping$SNP, pval=for_clumping$pval.exposure),
                           clump_kb=(window/1000),
                           clump_r2=clump_r2,
                           plink_bin=genetics.binaRies::get_plink_binary(),
                           bfile=opt$plink_binaries))
  if(inherits(out_clump, "try-error")){
    out_clump = try(TwoSampleMR::clump_data(dat = for_clumping, clump_kb = (window/1000), clump_r2 = clump_r2, pop = "EUR"))
    if(inherits(out_clump, "try-error")){
      stop("Problem with clumping!")
    }
    exposure_dat = exposure_format[which(exposure_format$SNP %in% out_clump$SNP),]
  } else {
    exposure_dat = exposure_format[which(exposure_format$SNP %in% out_clump$rsid),]
  }
#} else {
#  exposure_dat = exposure_format
#}
SNPS <- as.data.frame(exposure_format["SNP"])

message("\n\tClumping done!\n")
############################################################################## #
############################################################################## #
#### 4.3. Add F-statistic for the instrument/exposure ####
exposure_dat$samplesize.exposure<-opt$exp_ss_total
exposure_dat$r2 <- (exposure_dat$beta.exposure^2)/((exposure_dat$beta.exposure^2)+(exposure_dat$samplesize.exposure*(exposure_dat$se.exposure^2)))

exposure_dat$Fstat <- ((exposure_dat$samplesize.exposure[1] - length(exposure_dat$SNP) - 1)/length(exposure_dat$SNP))*(sum(exposure_dat$r2)/(1 - sum(exposure_dat$r2)))
#Columns to keep: "SNP";"effect_allele.exposure";"other_allele.exposure";"beta.exposure";"eaf.exposure";"pval.exposure";"se.exposure";"exposure";"r2";"Fstat"
#exposure_dat<-exposure_dat[,c("SNP","effect_allele.exposure","other_allele.exposure","beta.exposure","eaf.exposure","pval.exposure","se.exposure","exposure","r2","Fstat")]
fwrite(exposure_dat, paste0("/mnt/sda/gobemi01/PheWAS/PheWAS_MR_ANGPTL3/data/Fstat_PheWAS_Instruments_",opt$exp_pheno,"_",suff,"_pval.", pval_threshold,"_w",opt$window, "_r2.",opt$r2,".csv"), row.names = FALSE )
}
############################################################################## #
############################################################################## #
#### 5. Preparing results ####
{
  mr_res <- data.frame(matrix(nrow = 0, ncol = 14))
  colnames(mr_res) <- c("outcome", "exposure", "method", "top_snp", 
                        "b", "b_CIlow", "b_CIhigh", "se", "pval",
                        "r2_exp", "r2_out", "correct_causal_direction", 
                        "steiger_test_pval", "F_stat_b2_se2")
}
save_res <- FALSE

############################################################################## #
############################################################################## #
#### 6. Begin loop on each phenotype ####
{
  start_all <- Sys.time()
  if(use_from_to  & (from_na+to_na==2)){
    loop_command <- paste0(opt$from, ":", opt$to)
    if(opt$to >= nrow(phewas_dataset)){
      opt$to <- nrow(phewas_dataset)
    }
  } else {
    loop_command <- "1:nrow(phewas_dataset)"
  }
  for(i in eval(parse(text = loop_command))){
    message(paste0("\n\tProcessing phenotype ", opt$study, " (", phewas_dataset$outcome_name[i], ") (num :", i,") ...\n"))
    start <- Sys.time()
    if((i %% opt$save_every) == 0){
      save_res <- TRUE
    }
    opt$out_snp_col <- phewas_dataset$rsid_column[i]
    opt$out_chr_col <- phewas_dataset$chr_column[i]
    opt$out_pos_col <- phewas_dataset$pos_column[i]
    opt$out_ea_col <- phewas_dataset$effect_allele_column[i]
    opt$out_oa_col <- phewas_dataset$other_allele_column[i]
    opt$out_beta_col <- phewas_dataset$beta_column[i]
    opt$out_se_col <- phewas_dataset$se_column[i]
    opt$out_eaf_col <- phewas_dataset$eaf_column[i]
    opt$out_pval_col <- phewas_dataset$pvalue_column[i]
    opt$out_ss_col <- phewas_dataset$ss_column[i]
    opt$out_ss_total <- phewas_dataset$total_samplesize[i]
    opt$out_chr_col <- phewas_dataset$chr_column[i]

    window <- as.numeric(opt$window)
    
    if(opt$check_memory){
      test <- as.data.frame(system2("free", " --mega", stdout = TRUE))
      mem.col <- strsplit(test[1,], " +")[[1]]
      mem.total <- as.numeric(strsplit(test[2,], " +")[[1]][which(mem.col=='total')])
      mem.curr <- as.numeric(strsplit(test[2,], " +")[[1]][which(mem.col=='used')])
      mem.used <- as.numeric(sapply(strsplit(as.character(pryr::mem_used()), ' '), `[[`, 1))/1e06
      if(((mem.curr / mem.total) >= opt$check_memory_lim)){
        message(paste0("Warning : memory usage above ", opt$check_memory_lim*100 ," % . This script represents ", (mem.used/mem.curr)*100, "% of current memory usage."))
        if(((mem.curr / mem.total) >= 0.95) | ((mem.used/mem.curr) > opt$check_memory_lim)){
          stop("Memory usage too high. Stopping analysis.")
          gc()
          quit(save = "no", status=1)
        }
      } else {
        message(paste0("\n\t\t Script/Current memory usage : \t", (mem.used/mem.curr)*100, " %;\n\t\t Current/Total memory usage : \t", (mem.curr/mem.total)*100, " % "))
      }
    }
############################################################################## #
############################################################################## #
#### 6.1. Reading and preparing outcome file ####
{
  out_gwas <- as.data.frame(fread(phewas_dataset$path_GWAS[i]))
  #colnames(out_gwas) <- c(opt$out_chr_col,opt$out_pos_col,opt$out_ea_col,opt$out_oa_col,opt$out_snp_col,opt$out_pval_col,opt$out_beta_col,opt$out_se_col,opt$out_eaf_col)
  out_gwas$phenotype <- phewas_dataset$outcome[i]
  if(!is.na(phewas_dataset$total_samplesize[i])){
  out_gwas$samplesize <- as.numeric(phewas_dataset$total_samplesize[i])
  }else{
  out_gwas$samplesize <- NA
  }
  command_format <- paste0("outcome_format <- TwoSampleMR::format_data(dat = out_gwas,snps=exposure_format$SNP ,type = 'outcome'",
                           ", snp_col = '", opt$out_snp_col,
                           "', beta_col = '", opt$out_beta_col,
                           "', se_col = '", opt$out_se_col,
                           "', effect_allele_col = '", opt$out_ea_col,
                           "', other_allele_col = '", opt$out_oa_col,
                           "', pval_col = '", opt$out_pval_col,
                           "', phenotype_col = 'phenotype'")
  if(!is.na(opt$out_pos_col)){
    command_format <- paste0(command_format, ", pos_col = '", opt$out_pos_col, "'")
  }
  if(!is.na(opt$out_chr_col)){
    command_format <- paste0(command_format, ", chr_col = '", opt$out_chr_col, "'")
  }
  if(!is.na(opt$out_eaf_col)){
    command_format <- paste0(command_format, ", eaf_col = '", opt$out_eaf_col, "'")
  }
  if(!is.na(opt$out_ss_col)){
    command_format <- paste0(command_format, ", samplesize_col = '", opt$out_ss_col, "'")
  }
  command_format <- paste0(command_format, ")")
  eval(parse(text = command_format))
  # if(is.na(opt$out_gwas$samplesize) & !is.na(phewas_dataset$ncase[i])){
  #   outcome_format$samplesize.outcome <- (as.numeric(phewas_dataset$ncase[i]) + as.numeric(phewas_dataset$ncontrol[i]))
  # }else{
  #   outcome_format$samplesize.outcome <- NA  
  # }
  if(opt$out_p_is_log){
    outcome_format$pval.outcome <- 10**(outcome_format$pval.outcome)
  }
}
############################################################################## #
############################################################################## #
#### 6.3. Harmonization ####
message("\n\t Performing harmonization...\n")
dat_all <- TwoSampleMR::harmonise_data(exposure_dat = exposure_dat,
                                       outcome_dat = outcome_format)

#### 6.3.1. Saving MR dat (optional) ####
if(opt$save_dat){
  fwrite(x = dat_all, file = paste0(outname_dat, opt$exp_pheno,"_",opt$cis_trans,"_",phewas_dataset$outcome[i], "_SNPs_p.", pval_threshold, "_r2.", clump_r2,".txt"), sep = "\t")
}
############################################################################## #
############################################################################## #
#### 6.2. Proxies ####
if((nrow(exposure_dat) != nrow(dat_all)) & (opt$proxies)){
  message("Looking for proxies...")
  snps_exp_not_out <- subset(x = exposure_dat, subset = !(SNP %in% dat_all$SNP), select = "SNP")[[1]]
  success_proxies <- FALSE
  ## If you want to use files from server (remote), set path = ""
  gwasvcf::set_plink(path = "/usr/local/bin/plink")
  if(gwasvcf::check_plink()){
    proxies <- try(gwasvcf::get_ld_proxies(rsid = snps_exp_not_out, bfile = opt$plink_binaries, tag_r2 = 0.8, threads = opt$nthreads, out = paste0("./temp_files/", opt$study,"/",opt$exp_pheno, "/proxies.temp.", outcome_format$outcome[1], ".txt")))
    if(inherits(proxies, "try-error")){
      message("No proxies left to be found in reference file (local). No proxies will be added.")
      gwasvcf::set_plink(path = "")
      if(gwasvcf::check_plink()){
        proxies <- try(gwasvcf::get_ld_proxies(rsid = snps_exp_not_out, bfile = opt$plink_binaries, tag_r2 = 0.8, threads = opt$nthreads, out = paste0("./temp_files/", opt$study,"/",opt$exp_pheno, "/proxies.temp.", outcome_format$outcome[1],".txt")))
        if(inherits(proxies, "try-error")){
          message("No proxies left to be found in reference file (remote). No proxies will be added.")
        } else {
          success_proxies <- TRUE
        }
      }
    } else {
      success_proxies <- TRUE
    }
  } else {
    message("Finding proxies : gwasvcf::check_plink() is not true. Make sure you specify the plink software to use with gwasvcf::set_plink(). No proxies will be added.")
  }
  
  if(success_proxies){
    proxies$R2 <- proxies$R**2
    proxies <- proxies[with(proxies, order(R**2, decreasing = TRUE)),]
    # Trying top proxies
    if(use_from_to & (from_na+to_na==2)){
    proxies.top <- subset(x = proxies, subset = !duplicated(SNP_A))
    proxies.top <- subset(x = proxies.top, subset = SNP_B %in% outcome_format$SNP)
    if(nrow(proxies.top) == 0){
      proxies.top <- subset(x = proxies, subset = SNP_B %in% outcome_format$SNP)
      proxies.top <- proxies.top[with(proxies.top, order(R**2, decreasing = TRUE)),]
      proxies.top <- subset(x = proxies.top, subset = !duplicated(SNP_A))
    }
    }else{
      proxies.top <- subset(x = proxies, subset = !duplicated(SNP_A))
      proxies.top <- subset(x = proxies.top, subset = SNP_B %in% outcome_format$SNP)
      if(nrow(proxies.top) == 0){
        proxies.top <- subset(x = proxies, subset = SNP_B %in% outcome_format$SNP)
        proxies.top <- proxies.top[with(proxies.top, order(R**2, decreasing = TRUE)),]
        proxies.top <- subset(x = proxies.top, subset = !duplicated(SNP_A))  
      }
    }
    for(j in nrow(proxies.top)){
      exposure_dat[which(exposure_dat$SNP == proxies.top$SNP_A[j]), 
                   c("SNP", 
                     "effect_allele.exposure", 
                     "other_allele.exposure", 
                     "eaf.exposure", 
                     "beta.exposure", 
                     "pos.exposure")] <- c(
                       proxies.top$SNP_B[j], 
                       proxies.top$B1[j], 
                       proxies.top$B2[j], 
                       proxies.top$MAF_B[j], 
                       ifelse(test = proxies.top$R[j] < 0, 
                              yes = -1*exposure_dat$beta.exposure, 
                              no = exposure_dat$beta.exposure), 
                       proxies.top$BP_B[j])
    }
    dat_all <- TwoSampleMR::harmonise_data(exposure_dat = exposure_dat,
                                           outcome_dat = outcome_format)
    if(opt$save_dat){
      fwrite(x = dat_all, file = paste0(outname_dat,opt$exp_pheno,"_",opt$cis_trans,"_",outcome_format$outcome[1] , "_SNPs_p.", pval_threshold, "_r2.", clump_r2,".txt"), sep = "\t")
    }
  }
}
if(file.exists(paste0("./temp_files/", opt$study,"/",opt$exp_pheno, "/proxies.temp.", outcome_format$outcome[1],".txt"))){
  unlink(paste0("./temp_files/", opt$study,"/",opt$exp_pheno, "/proxies.temp.", outcome_format$outcome[1], ".txt"))
}
if(!opt$reverse){
  dat_all$pos.outcome <- dat_all$pos.exposure
}else {
  dat_all$pos.exposure <- dat_all$pos.outcome
}
if(nrow(dat_all)==0) {
  message("\t\t *** WARNING : no SNPs in the european reference panel. \n\t\t NO HARMONIZATION WILL BE MADE ***")
  message("No SNP found in reference panel")
  next
}
if(opt$singleSNP_multiSNP == "singleSNP"){
  min_pval <- min(dat_all$pval.exposure)
  dat_all <- dat_all %>% filter(dat_all$pval.exposure == min_pval)
}  
dat <- dat_all
maf <- ifelse(dat$eaf.exposure<=0.5,dat$eaf.exposure,1-dat$eaf.exposure)
############################################################################## #
############################################################################## #
#### 6.4. Perform MR ####
if(nrow(dat) == 1){
  mr_wald <- mr_wald_ratio(b_exp = dat$beta.exposure, b_out = dat$beta.outcome, se_exp = dat$se.exposure, se_out = dat$se.outcome)
  if(length(mr_wald) == 0){
    message("No SNPs available for MR for ", outcome_format$outcome[1], "\n")
    unprocessed_proteins <- rbind(unprocessed_proteins, data.frame(prot = outcome_format$outcome[1], reason = as.character("MR : No SNPs remaining for MR")))
    fwrite(x = unprocessed_proteins,
           file = unprocessed_proteins_file,
           sep = "\t", append = FALSE)
    next
  }
  CIs <- generate_odds_ratios(mr_wald)
  mr_prot <- data.frame(outcome=dat$outcome[1], exposure=dat$exposure[1], method="Wald ratio", top_snp = dat$SNP, b = mr_wald$b, b_CIlow = CIs$lo_ci, b_CIhigh = CIs$up_ci, se = mr_wald$se, pval = mr_wald$pval)
  #rm(mr_wald)
  
}
if(nrow(dat) == 0){
  message("\n\tNo SNPs remaining after MR. Skipping protein...\n")
  next
}

############################################################################## #
############################################################################## #
#### 7. Perform additional analyzes ####

## 7.1. MR Steiger directionality test ####
if(nrow(dat) >= 1){
  dat$r_exp <- ifelse(test = !is.na(opt$exp_ncase),
                      yes = get_r_from_lor(lor=dat$beta.exposure,af=dat$eaf.exposure,ncase=as.numeric(opt$exp_ncase),
                                           ncontrol=(as.numeric(opt$exp_ss_total)-as.numeric(opt$exp_ncase)),prevalence=(as.numeric(opt$exp_ncase) / as.numeric(opt$exp_ss_total))),
                      no = get_r_from_pn(dat$pval.exposure, as.numeric(opt$exp_ss_total)))	# for Steiger
  
  dat$r_out <- ifelse(test = !is.na(phewas_dataset$ncase[i]),
                      yes = get_r_from_lor(lor=dat$beta.outcome,af=maf,ncase=as.numeric(phewas_dataset$ncase[i]),
                                           ncontrol=as.numeric(phewas_dataset$ncontrol[i]),prevalence=as.numeric(phewas_dataset$ratio[i])),
                      no = get_r_from_pn(dat$pval.outcome,as.numeric(opt$out_ss_total)))
  
  steiger <- mr_steiger(
    p_exp = dat$pval.exposure, 
    p_out = dat$pval.outcome, 
    n_exp = as.numeric(opt$exp_ss_total), 
    n_out = as.numeric(opt$out_ss_total), 
    r_xxo = 1, 
    r_yyo = 1,
    r_exp = dat$r_exp,
    r_out = dat$r_out
  )
  mr_prot$r2_exp <- steiger[c(1)]
  mr_prot$r2_out <- steiger[c(2)]
  mr_prot$correct_causal_direction <- steiger[c(5)]
  mr_prot$steiger_test_pval <- steiger[c(6)]
## 7.2. Adding F-stat ####
  F_stat_b2_se2 = function(beta, se) {sum(beta^2/se^2)}
  mr_prot$F_stat_b2_se2 = F_stat_b2_se2(beta=dat$beta.exposure,se=dat$se.exposure)
}else{
  mr_prot$r2_exp <- NA
  mr_prot$r2_out <- NA
  mr_prot$correct_causal_direction <- NA
  mr_prot$steiger_test_pval <- NA
  mr_prot$F_stat_b2_se2 = NA 
}
############################################################################## #
############################################################################## #
#### 8. Save results ####
if(length(mr_prot)==9){
mr_res <- rbind(mr_res, mr_prot)
} else {
mr_prot <- mr_prot[3:14]  
mr_res <- rbind(mr_res, mr_prot)
}
if(((i %% opt$save_every) == 0) | save_res){
  fwrite(x = mr_res,
         file = outname,
         sep = "\t", append = FALSE)
  save_res <- FALSE
}
# rm(list = c("pwas_gene", "exposure", "mr_prot", "dat", "exposure_format"))
end <- Sys.time()
print(end - start)
}
fwrite(x = mr_res,
       file = outname,
       sep = "\t", append = FALSE)
end_all <- Sys.time()
message("Total time to run loop on all specified phenotypes : \n")
print(end_all - start_all)
#rm(list=ls())
message("==================== DONE ==================== \n")

}
unlink(paste0("./temp_files/",opt$study,"/",opt$exp_pheno,"/",opt$study,"_snps_list_from.", opt$from, "_to.", opt$to, ".txt" ))
unlink(paste0("./temp_files/",opt$study,"/",opt$exp_pheno,"/",opt$study,"_gwas_from.", opt$from, "_to.", opt$to, ".txt"))

