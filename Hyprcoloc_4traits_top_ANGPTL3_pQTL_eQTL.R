#!/usr/bin/env Rscript

#Terminal command
#for top eQTL:
#nohup ./hyprcoloc_plot_4traits_top_ANGPTL3_pQTL_eQTL.R --gwas_table=/mnt/sda/gobemi01/Data_instruments/gwas_table_cardiometabolic_06.03.2023.xlsx --chr 1 --gene ANGPTL3 --pQTL_eQTL eQTL --hyprcoloc=TRUE --out /mnt/sda/gobemi01/ANGPTL3/HyPrColoc/top_eQTL_pQTL_ANGPTL3_TG_APOB_CAD/ --gene_ensembl ENSG00000132855.5 --top_snp rs11207978 --top_pos 62514841 >./nohupEQTL.txt &
#for top pQTL:
#nohup ./hyprcoloc_plot_4traits_top_ANGPTL3_pQTL_eQTL.R --gwas_table=/mnt/sda/gobemi01/Data_instruments/gwas_table_cardiometabolic_06.03.2023.xlsx --chr 1 --gene ANGPTL3 --pQTL_eQTL pQTL --hyprcoloc=TRUE --out /mnt/sda/gobemi01/ANGPTL3/HyPrColoc/top_eQTL_pQTL_ANGPTL3_TG_APOB_CAD/ --gene_ensembl ENSG00000132855.5 --top_snp rs10889352 --top_pos 62633352 >./nohupPQTL.txt &

# Or run it in R:
#for top eQTL:
opt<- data.frame(gwas_table="/mnt/sda/gobemi01/Data_instruments/gwas_table_cardiometabolic_06.03.2023.xlsx",chr=1,gene="ANGPTL3",pQTL_eQTL="eQTL", hyprcoloc=TRUE,out="/mnt/sda/gobemi01/ANGPTL3/HyPrColoc/top_eQTL_pQTL_ANGPTL3_TG_APOB_CAD/",gene_ensembl = "ENSG00000132855.5",top_snp="rs11207978",top_pos=62514841)
#for top pQTL:
#opt<- data.frame(gwas_table="/mnt/sda/gobemi01/Data_instruments/gwas_table_cardiometabolic_06.03.2023.xlsx",chr=1,gene="ANGPTL3",pQTL_eQTL="pQTL" ,hyprcoloc=TRUE,out="/mnt/sda/gobemi01/ANGPTL3/HyPrColoc/top_eQTL_pQTL_ANGPTL3_TG_APOB_CAD/",gene_ensembl = "ENSG00000132855.5",top_snp="rs10889352",top_pos=62633352)

library(data.table)
library(ggrepel)
library(grid)
library(gridExtra)
library(gtable)
library(R.devices)
library(tidyverse)
library(optparse)
source("/mnt/sda/boujer01/Scripts/my_stack_assoc_plot.R")

#### Function ####
library(grid)
mygenes <- readRDS(file = "/mnt/sda/boujer01/Scripts/gassoc_myGenes")
source("/mnt/sda/gobemi01/Functions/my_stack_assoc_plot_JéJé_copy.R")
source("/mnt/sda/gobemi01/Gobeil.R/R/gassocplot_colorcustom.R")
option_list = list(
  make_option("--chr", action="store", default=1, type='numeric',
              help=" Chromosome to analyze "),
  make_option("--gwas_table", action="store", default=NA, type='character',
              help=" gwas_table to loop "),
  make_option("--gene", action="store", default=NA, type='character',
              help=" Gene to analyze (must correspond to gene_id or gene.id in the MRdat file name) "),
  make_option("--gene_ensembl", action="store", default=NA, type='character',
              help=" Gene ENSEMBL to analyze "),
  make_option("--pQTL_eQTL", action="store", default=NA, type='character',
              help="Top pQTL or eQTL for Trait4?"),
  make_option("--top_snp", action="store", default=NA, type='character',
              help="Top SNP rsid to highlight in plot"),
  make_option("--top_pos", action="store", default=NA, type='numeric',
              help="Top SNP position b38 to highlight in plot"),
  make_option("--hyprcoloc", action="store", default=TRUE, type='logical',
              help="Perform HyPrColoc ? "),
  make_option("--out", action="store", default=NA, type='character',
              help=" Output directory (full path) for regional plots without group name and chromosomes folders (eg. : ./regional_plots/2_groups/) "))
opt = parse_args(OptionParser(option_list=option_list))

if(is.na(opt$gene)){stop("No gene was specified for the analysis.")}
if(is.na(opt$out)){stop("No output directory was specified for the analysis.")}

gene_name <- opt$gene
gene <- opt$gene
gene_chr <- opt$chr
chr <- opt$chr
out <- opt$out
pQTL_eQTL <- opt$pQTL_eQTL
window=1000000
n_traits <- 4
##


wd = paste0(opt$out)
if(!dir.exists(wd)){
  dir.create(wd)
}
setwd(wd)

gene.region <- fread("/mnt/sda/girarn01/vte/proteins_decode.txt") 
gene.region <- gene.region %>% filter(prot_name == "ANGPTL3")
gene.region$chr_decode <- paste0("chr", gene.region$CHR)

## TRAIT 1 : CAD
trait1 <-as.data.frame(fread("/mnt/sda/gobemi01/GWAS/GWAS_in_b38/b38_GCST90132314_buildGRCh37_rsids.tsv",header = TRUE, stringsAsFactors = FALSE))
trait1 <- subset(x = trait1, select = c("rsid", "beta", "standard_error", "effect_allele", "other_allele", "chromosome", "pos_b38", "p_value","effect_allele_frequency"))
colnames(trait1) <- c("SNP", "beta", "se", "ea", "oa", "chr", "pos", "pval","eaf")
trait1 <- trait1 %>% filter(chr == gene.region$CHR & between(pos, gene.region$start-window/2, gene.region$end+window/2) == T)
trait1$pheno <- paste0("CAD")
trait1$SS <- 1165720

## TRAIT 2 : ApoB ####
trait2 <-fread('/mnt/sda/gobemi01/GWAS/GWAS_in_b38/b38_apob_30640_mod.csv.gz')
trait2 <- subset(x = trait2, select = c("rsid", "beta", "se", "alt_allele", "ref_allele", "chr", "pos_b38", "pval","EAF"))
colnames(trait2) <- c("SNP", "beta", "se", "ea", "oa", "chr", "pos", "pval","eaf")
trait2 <- trait2 %>% filter(chr == chr & between(pos, gene.region$start-window/2, gene.region$end+window/2) == T)
trait2 <- trait2 %>% filter(SNP %in% trait1$SNP) %>% .[order(.$pval),] %>% filter(!duplicated(SNP)) 
trait2$pheno <- paste0("ApoB")
trait2$SS <- 439214

## TRAIT 3 : TG ####
trait3 <- fread("/mnt/sda/gobemi01/GWAS/GWAS_in_b38/b38_logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz")
trait3 <- subset(x = trait3, select = c("rsID", "EFFECT_SIZE", "SE", "ALT", "REF", "CHROM", "pos_b38", "pvalue","POOLED_ALT_AF"))
colnames(trait3) <- c("SNP", "beta", "se", "ea", "oa", "chr", "pos", "pval.txt","eaf")
trait3 <- trait3 %>% filter(chr == chr & between(pos, gene.region$start-window/2, gene.region$end+window/2) == T)
# Given the very small pval values that reach a floor, we'll use the Rmpfr package to generate a S4 class.
library(Rmpfr)
trait3$pval <- mpfr(x=trait3$pval.txt,precBits=15)
trait3 <- trait3[which(trait3$SNP %in% trait1$SNP),]
trait3 <- trait3[BiocGenerics::order(trait3$pval),]
trait3 <- trait3[which(!duplicated(trait3$SNP)),]
trait3$pheno <- paste0("TG")
trait3$SS <- 1320016

if(opt$pQTL_eQTL=="pQTL"){
  ### TRAIT 4 : ANGPTL3 (pQTL) ####
  trait4 <- as.data.frame(fread(paste0("/mnt/sda/dbs_Web/deCODE/10391_1_ANGPTL3_ANGL3.txt.gz"), header = TRUE, stringsAsFactors = FALSE))
  trait4 <- subset(x = trait4, select = c("rsids", "Beta", "SE", "effectAllele", "otherAllele", "Chrom", "Pos", "Pval","ImpMAF"))
  trait4$Chrom <- gsub(x = trait4$Chrom, pattern = 'chr', replacement = '')
  trait4$Chrom <- as.numeric(trait4$Chrom)
  colnames(trait4) <- c("SNP", "beta", "se", "ea", "oa", "chr", "pos", "pval","eaf")
  trait4 <- trait4 %>% filter(chr == chr & between(pos, gene.region$start-window/2, gene.region$end+window/2) == T)
  trait4 <- trait4 %>% na.omit(.) %>% filter(SNP %in% trait1$SNP) %>% .[order(.$pval),] %>% filter(!duplicated(SNP))
  trait4$pheno <- paste0("ANGPTL3.pQTL")
  trait4$SS <- 35559
} 
if(opt$pQTL_eQTL=="eQTL"){
  ### TRAIT 4 : ANGPTL3 (eQTL) ####
  gencode = as.data.frame(data.table::fread("/mnt/sda/couchr02/Eloi/tensorqtl/gencode.v34.GRCh38.genes.txt"))
  gencode_small = gencode[which(gencode$gene_name=="ANGPTL3"),]
  gencode_small$chr = as.numeric(gencode_small$chr)
  gene_id = unique(gencode_small$gene_id)
  chrom = unique(gencode_small$chr)
  mystart = min(gencode_small$start)
  myend = max(gencode_small$end)
  setwd("/mnt/sda/couchr02/Eloi/tensorqtl/")
  trait4 <- as.data.frame(fread(cmd=paste0("tabix -h /mnt/sda/couchr02/Eloi/tensorqtl/tensorqtl_ref_all.cis_qtl_pairs.genes.rsids.ALL_ASSOCIATIONS_sorted.txt.gz ",chrom,":",mystart-window/2,"-",myend+window/2," | grep \"phenotype_id\\|",gene_id,"\"")))
  colnames(trait4) <- c("chr","pos","SNP","gene","ensembl_id","eaf","pval","beta","se","ea","oa")
  trait4$pval <- as.numeric(trait4$pval)
  trait4 <- trait4 %>% filter(SNP %in% trait1$SNP) %>% .[order(.$pval),] %>% filter(!duplicated(SNP))
  trait4$pheno <- paste0("ANGPTL3.eQTL")
  trait4$SS <- 246
}
setwd(wd)
traits <- c(trait1$pheno[1],trait2$pheno[1],trait3$pheno[1],trait4$pheno[1])

words1 <- "awk"
args1 <- paste0("'$2==", gene_chr, " {print}' /mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b38_rsid_maf_small.txt >> ", out, "snps_ref_chr", gene_chr, ".txt")
system2(words1, args1)

snps_ref <- as.data.frame(fread(file = paste0(out, "snps_ref_chr", gene_chr, ".txt"), header = FALSE, stringsAsFactors = FALSE))
colnames(snps_ref) <- c("rsid", "chr", "pos_b38", "a0", "a1", "EUR", "eaf")
snps_ref <- subset(x = snps_ref, select = c("rsid", "chr", "pos_b38"))
file.remove(paste0(out, "/snps_ref_chr", gene_chr, ".txt"))


for(i in 1:n_traits){
  eval(parse(text = paste0("trait", i, " <- subset(x = trait", i, ", subset = chr == opt$chr)")))
}

if(!file.exists(paste0("chr",opt$chr,".bed")) | !file.exists(paste0("chr",opt$chr,".bim")) | !file.exists(paste0("chr",opt$chr,".fam")) | !file.exists(paste0("chr",opt$chr,".nosex"))){
  system2(
    "plink",
    args = c(
      "--vcf",
      paste0("/home/couchr02/Mendel_Commun/Nicolas/1000G_phase3/ALL.chr",opt$chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"),
      "--chr",
      opt$chr,
      "--keep",
      "/home/couchr02/Mendel_Commun/Nicolas/1000G_phase3/List_EUR_non_related.txt",
      "--make-bed",
      "--out",
      paste0("chr",opt$chr)
    )
  )
}

gene <- opt$gene
message(paste0("\n *** Processing gene : ", gene, "*** \n"))

### Creating LD matrix
top_eqtl <- opt$top_snp
top_pos<-opt$top_pos
top_chr <- opt$chr




for(i in 1:n_traits){
  if(i==3){
    eval(parse(text = paste0("trait", i, " <- trait", i, "[BiocGenerics::order(trait", i, "$pval, decreasing = FALSE), ]")))
    eval(parse(text = paste0("trait", i, " <- BiocGenerics::subset(x = trait", i, ", subset = !duplicated(SNP))")))
    eval(parse(text = paste0("trait", i, " <- BiocGenerics::subset(x = trait", i, ", subset = (!is.na(beta)) & (!is.na(se)))")))
    }else{
  eval(parse(text = paste0("trait", i, " <- trait", i, "[order(trait", i, "$pval, decreasing = FALSE), ]")))
  eval(parse(text = paste0("trait", i, " <- subset(x = trait", i, ", subset = !duplicated(SNP))")))
  eval(parse(text = paste0("trait", i, " <- subset(x = trait", i, ", subset = (!is.na(beta)) & (!is.na(se)))")))
} }

trait4 <-subset(x = trait4, subset = (SNP %in% trait1$SNP) & (SNP %in% trait2$SNP) & (SNP %in% trait3$SNP)) #& (SNP %in% trait5$SNP))
trait1 <-subset(x = trait1, subset = SNP %in% trait4$SNP)
trait2 <-subset(x = trait2, subset = SNP %in% trait4$SNP)
trait3 <-BiocGenerics::subset(x = trait3, subset = SNP %in% trait4$SNP)


if(opt$hyprcoloc){
  for(i in 1:n_traits){
    type <- ifelse(test = i == 4, yes = 'exposure', no = 'outcome')
    if(i==3){
    form3 <- trait3
    colnames(form3) <- c("SNP","beta.outcome","se.outcome","effect_allele.outcome","other_allele.outcome","chr.outcome","pos.outcome","pval.txt","eaf.outcome","pval.outcome","outcome","samplesize.outcome")
    form3$mr_keep.outcome <- TRUE
    form3$pval_origin.outcome <- "reported"
    form3$id.outcome <- "Trait3"
  }else{
    eval(parse(text = paste0("form", i, " <- TwoSampleMR::format_data(dat = trait", i, ", type = '", type, "',min_pval = 1e-600 ,snp_col = 'SNP', beta_col = 'beta', se_col = 'se', effect_allele_col = 'ea', other_allele_col = 'oa', eaf_col='eaf' ,chr_col = 'chr', pos_col = 'pos', pval_col = 'pval', samplesize_col='SS',phenotype_col='pheno')" )))
    }
    }
  for(i in 1:3){
    eval(parse(text = paste0("harm4_", i, " <- TwoSampleMR::harmonise_data(exposure_dat = form4, outcome_dat = form", i, ")")))
    if(i == 3){
      eval(parse(text = paste0("harm4_3 <- subset(x = harm4_3, subset = (SNP %in% harm4_1$SNP) & (SNP %in% harm4_2$SNP))")))
      for(j in 1:2){
        eval(parse(text = paste0("harm4_", j, " <- subset(x = harm4_", j, ", subset = (SNP %in% harm4_3$SNP))")))
      }
    }
  }
  
  betas <- data.frame(matrix(nrow=nrow(harm4_3), ncol = n_traits))
  ses <- data.frame(matrix(nrow=nrow(harm4_3), ncol = n_traits))
  betas$X4 <- harm4_3$beta.exposure
  ses$X4 <- harm4_3$se.exposure
  for(i in 1:3){
    for(j in c("beta", "se")){
      eval(parse(text = paste0(j, "s$X", i, " <- harm4_", i, "$", j, ".outcome")))
      if(i == 3){
        eval(parse(text = paste0("colnames(", j, "s) <- traits; ", j, "s <- as.matrix(", j, "s); rownames(", j, "s) <- harm4_3$SNP")))
      }
    }
  }
  rsids <- rownames(betas)
  
  hypr_coloc_res <- hyprcoloc::hyprcoloc(betas, ses, trait.names=traits, snp.id=rsids, binary.outcomes = c(1,0,0,0),snpscores = TRUE) 
  
  #top_eqtl <- hypr_coloc_res[["results"]][["candidate_snp"]]
  
  png(filename = paste0("res_sens_plot_4traits_top_ANGPTL3_",pQTL_eQTL,"_CAD_ApoB_TG.png"), type = 'cairo')
  hypr_coloc_plot <- hyprcoloc::sensitivity.plot(betas, ses, trait.names = traits, similarity.matrix = TRUE)
  #hypr_coloc_plot
  dev.off()
  #hypr_coloc_corr <- hyprcoloc::cred.sets(hypr_coloc_res)
  fwrite(x = hypr_coloc_res$results, file = paste0("hyprcoloc_res_4traits_top_ANGPTL3_",pQTL_eQTL,"_CAD_ApoB_TG.txt"), sep = "\t")
  
}

fwrite(as.data.frame(trait4$SNP),"list_snps.txt",row.names=FALSE,col.names=FALSE)

system2(
  "plink",
  args = c(
    "--bfile",
    paste0("chr",opt$chr),
    "--extract",
    "list_snps.txt",
    "--make-bed",
    "--out",
    "1000G_snp"
  )
)
system2(
  "plink",
  args = c(
    "--bfile",
    "1000G_snp",
    "--r square",
    "--out",
    "LD_Matrix_1000G"
  )
)
if(!file.exists("LD_Matrix_1000G.ld")){
  message("\n *** No SNPs could be extracted with PLINK. Skipping gene... *** \n")
  file.remove("LD_Matrix_1000G.log")
  file.remove("1000G_snp.log")
  file.remove("1000G_snp.nosex")
  file.remove("list_snps.txt")
  next
}
LD_Mat = read.table("LD_Matrix_1000G.ld", header = FALSE, stringsAsFactors = F)
bim_file = read.table("1000G_snp.bim", header = FALSE, stringsAsFactors = F)
if(sum(duplicated(bim_file$V2)) > 0){
  to_remove <- which(duplicated(bim_file$V2))
  bim_file <- bim_file[-c(to_remove),]
  LD_Mat <- LD_Mat[-c(to_remove),-c(to_remove)]
}
colnames(LD_Mat) = bim_file$V2
rownames(LD_Mat) = bim_file$V2

if(nrow(bim_file) > 1){
  bim_file <- bim_file[!(duplicated(bim_file$V2)),]
  LD_Mat <- LD_Mat[!(duplicated(colnames(LD_Mat))), !(duplicated(colnames(LD_Mat)))]
} else {
  message("\n *** Not enough SNPs kept for LD Matrix (nsnps < 2). Skipping... *** \n")
  file.remove("LD_Matrix_1000G.log")
  file.remove("1000G_snp.log")
  file.remove("1000G_snp.nosex")
  file.remove("list_snps.txt")
  next
}

file.remove("1000G_snp.log")
file.remove("1000G_snp.nosex")
file.remove("1000G_snp.bed")
file.remove("1000G_snp.bim")
file.remove("1000G_snp.fam")

file.remove("LD_Matrix_1000G.ld")
file.remove("LD_Matrix_1000G.log")
file.remove("LD_Matrix_1000G.nosex")

file.remove("list_snps.txt")

### Creating plots
for(i in 1:n_traits){
  eval(parse(text = paste0("trait", i, " <- trait", i, "[order(trait", i, "$SNP), ]")))
}
z_mat <- data.frame(matrix(nrow=nrow(trait4), ncol = n_traits))
markers_mat <- data.frame(marker = trait4$SNP, chr = trait4$chr, pos = trait4$pos)

for(i in 1:n_traits){
  eval(parse(text = paste0("trait", i, "$pval <- as.numeric(trait", i, "$pval)")))
  eval(parse(text = paste0("z_mat[,",i, "] <- abs(trait", i, "$beta/trait", i,"$se)")))
}
rownames(z_mat) = markers_mat$marker
markers_mat <- markers_mat[which(markers_mat$marker %in% colnames(LD_Mat)),]
assign("ldmat", LD_Mat[markers_mat$marker,markers_mat$marker])
markers_mat$pos <- as.integer(markers_mat$pos)
markers_mat$chr <- as.integer(markers_mat$chr)
z_mat <- subset(x = z_mat, subset = rownames(z_mat) %in% markers_mat$marker)

png(filename = paste0("/mnt/sda/gobemi01/ANGPTL3/HyPrColoc/top_eQTL_pQTL_ANGPTL3_TG_APOB_CAD/hyprcoloc_plot_4traits_top_ANGPTL3_",pQTL_eQTL,"_CAD_ApoB_TG.png"), width = 800, height = 1200, type = 'cairo')
mystack <- my_stack_assoc_plot(markers_mat,z_mat,ldmat, traits=traits, 
                               top.marker=top_eqtl, 
                               x.min = as.integer(gene.region$start)-(5e05), 
                               x.max = as.integer(gene.region$end)+(5e05),build = 38)
mystack

dev.off()

file.remove(paste0("chr",opt$chr,".log"))
file.remove(paste0("chr",opt$chr,".nosex"))
file.remove(paste0("chr",opt$chr,".bed"))
file.remove(paste0("chr",opt$chr,".bim"))
file.remove(paste0("chr",opt$chr,".fam"))

rm(list=ls())