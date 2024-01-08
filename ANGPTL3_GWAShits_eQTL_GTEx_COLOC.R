#!/usr/bin/env Rscript

wd_data = "/mnt/sda/gobemi01/ANGPTL3/coloc_res"
out = "/mnt/sda/gobemi01/ANGPTL3/coloc_res"
setwd(wd_data)

library(tidyverse)
library(data.table)
library(stringr)
library(TwoSampleMR)
library(coloc)
library(locuscomparer)
library(ggplot2)

genes = c("ANGPTL3","LPL","APOA5","APOA4","APOC3","APOA1","ANGPTL4","APOE") #c("ANGPTL3")  #c("ANGPTL3","LPL","APOA5-APOA4-APOC3-APOA1","ANGPTL4","APOE") 	# can be a list of genes
window = 1000000
locuscompare_plot = "Yes"  # "Yes" or "No"
target_snp = NULL #"rs773045547"	# a unique SNP with "rsname" or NULL
block_bp <- 10000000
name_gwas = "ANGPTL3_deCODE"


# selection of GWAS top hits

gwas <- as.data.frame(fread(paste0("/mnt/sda/dbs_Web/deCODE/10391_1_ANGPTL3_ANGL3.txt.gz"), header = TRUE, stringsAsFactors = FALSE))
gwas <- subset(x = gwas, select = c("rsids", "Beta", "SE", "effectAllele", "otherAllele","ImpMAF", "Chrom", "Pos", "Pval"))
gwas$Chrom <- gsub(x = gwas$Chrom, pattern = 'chr', replacement = '')
gwas$Chrom <- as.numeric(gwas$Chrom)
colnames(gwas) <- c("SNP", "beta", "se", "ea", "oa","maf", "chr", "pos", "pval")
gwas$pval <- as.numeric(gwas$pval)

# Loading annotation file with locus
mySNPs = as.data.frame(fread(paste0("/mnt/sda/gobemi01/Data_instruments/for_manhattan_plot/annot_ANGPTL3_deCODE_for_Manhattan.txt")))
mySNPs = mySNPs[-c(2),]
mySNPs$Locus = c("ANGPTL3","LPL","APOA5-APOA4-APOC3-APOA1","ANGPTL4","APOE")

# selection of tissues of interest from a list 

tissue_dir = list.files("/mnt/sda/GTEx_v8_Mendel_Commun/eQTL/GTEx_Analysis_v8_eQTL_expression_matrices/")
tissue_loop = sub(".v8.normalized_expression.bed.gz","",tissue_dir)
tissue_loop = tissue_loop[c(1,3,61,97)]  # Adipose_Subcutaneous, Adipose_Visceral_Omentum, Liver

# loading snps file containing rsid, chr and pos (build 38) 

snps = as.data.frame(fread("/home/couchr02/workspace/1000G_Phase3/1000G_Phase3_b38_rsid_maf_small.txt"))

# loading genes file in order to retreive ENSEMBL name and also starting and ending positions of the gene (build 38)

gencode = as.data.frame(fread("/home/couchr02/Mendel_Commun/Christian/GTEx_v8/gencode.v26.GRCh38.genes.txt"))

  # loading and format ANGPTL3 data

  results_tot = data.frame()
  for (gene in genes) {
    
    # retrieving some useful infos from genes file
    
    gencode_small = gencode[which(gencode$gene_name==gene),]
    gencode_small$chr = as.numeric(gencode_small$chr)
    gene_id = unique(gencode_small$gene_id)
    CHR = unique(gencode_small$chr)
    mystart = min(gencode_small$start)
    myend = max(gencode_small$end)
    
    # loading the associations file created by QTLtools and retrieving the rsids for the gene (by chr and pos)
    
    QTL_gene = as.data.frame(fread(paste0("/mnt/sda/gobemi01/ANGPTL3/coloc_res/QTLtools/All_tissues_",gene,"_QTLtools_nominal_1_window_1000kb.txt")))
    snps_gene = snps[which(snps$chr %in% QTL_gene$chr_top_variant & snps$pos_b38 %in% QTL_gene$start_top_variant),]
    
    # Exposure data and format
    
    exp_gwas <- gwas %>% filter(chr == CHR & between(pos, mystart-window/2, myend+window/2) == T)
    exp_gwas <- na.omit(exp_gwas) %>% unique(.) %>% .[order(.$pval),] %>% filter(!duplicated(SNP))
    exp_gwas$pheno <- paste0("ANGPTL3.pQTL")
    exp_gwas$SS <- 35559
    
    
    for (tissue in tissue_loop) {
      
      # selecting one tissue to analyze
      
      QTL = QTL_gene[which(QTL_gene$tissue==tissue),] 
      
      # merging rsids with gene (by chr and pos) and keeping the SNPs present in the gene and the outcome for the
      # exposure (eQTL)
      
      merged = merge(QTL,snps_gene,by.x=c("chr_top_variant","start_top_variant"),by.y=c("chr","pos_b38"),sort=F)
      
      # keeping the SNPs present in the gene and the outcome for the outcome
      s<-Gobeil.R::get_proxies(exp_gwas,merged,snp1_col="SNP",snp2_col="rsid",rsq = 0.8, query_file="/mnt/sda/gobemi01/queries/eqtl_combined_query_snp_list.txt")
     
    ### TRAIT exposure : ANGPTL3 (pQTL) ####
    exposure <- format_data(exp_gwas, 
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
                          min_pval = 1e-400
    ) %>% as.data.table(.) %>% .[,id.exposure := "ANGPTL3.pQTL"]
    
     # Outcome / GTEX tissue-specific eQTLs
      
      outcome <- format_data(dat = merged,
                              type = "outcome",
                              phenotype_col = "tissue",
                              snp_col = "rsid",
                              beta_col = "beta",
                              se_col = "se",
                              effect_allele_col = "a0",
                              other_allele_col = "a1",
                              pval_col = "p_val",
                             eaf_col = "EUR",
                              chr_col = "chr_top_variant",
                              pos_col = "start_top_variant",
                              gene_col = "probe",
                              samplesize_col = "sample_size",
                             min_pval = 1e-400
         ) %>% as.data.table(.) %>% .[,id.outcome := "GTEx.eQTL"]
      

      # harmonizing the two datasets to make the coloc analysis easier	
      
      coloc_dat = harmonise_data(exposure, outcome, action=2)
      if(nrow(coloc_dat)==0){
        next
      } else {
      
      ### COLOC ####
      
      coloc_dat$maf.exposure = ifelse(coloc_dat$eaf.exposure<=0.5,coloc_dat$eaf.exposure,1-coloc_dat$eaf.exposure)
      coloc_dat$maf.outcome = ifelse(coloc_dat$eaf.outcome<=0.5,coloc_dat$eaf.outcome,1-coloc_dat$eaf.outcome)
      
      dat_window = coloc_dat[which(coloc_dat$pos.exposure >= mystart-(window/2) & coloc_dat$pos.exposure <= myend+(window/2)),]
      
      samplesize_out = as.numeric(outcome$samplesize.outcome[1])
      
      dataset1 <- data.table(pvalues = dat_window$pval.exposure,
                             N = 35346,
                             MAF = dat_window$maf.exposure,
                             beta = dat_window$beta.exposure,
                             varbeta = (dat_window$se.exposure)^2,
                             type = "quant",
                             snp = dat_window$SNP) %>% na.omit(.)
      
      dataset2 <- data.table(pvalues = dat_window$pval.outcome,
                             N = samplesize_out,
                             MAF = dat_window$maf.outcome,
                             beta = dat_window$beta.outcome,
                             varbeta = (dat_window$se.outcome)^2,
                             type = "quant",
                             snp = dat_window$SNP) %>% na.omit(.)
      
      res_coloc <- coloc.abf(dataset1, dataset2)
      
      results = data.frame(unique(outcome$outcome),gene,tissue,res_coloc$summary[1],res_coloc$summary[2],res_coloc$summary[3],res_coloc$summary[4],res_coloc$summary[5],res_coloc$summary[6])
      colnames(results) = c("outcome","gene","tissue","nSNP","coloc_PP.H0","coloc_PP.H1","coloc_PP.H2","coloc_PP.H3","coloc_PP.H4")
      results_tot = bind_rows(results_tot,results)
      
      ### end of coloc ###
      
      ### locuscompare ###			
      
      if (locuscompare_plot=="Yes") {
        
        dat_exp = coloc_dat[,c("SNP","pval.exposure")]
        colnames(dat_exp)=c("rsid","pval")
        
        dat_out = coloc_dat[,c("SNP","pval.outcome")]					  
        colnames(dat_out)=c("rsid","pval")
        
        locus_zoom <- locuscompare(
          in_fn1 = dat_out,
          in_fn2 = dat_exp,
          genome = "hg38",
          snp = target_snp
        )
        
        ggsave(locus_zoom,file=paste0("/mnt/sda/gobemi01/ANGPTL3/coloc_res/",name_gwas,"_",gene,"_",target_snp,"_",tissue,".png"),width=20, height=15,units="cm",dpi=200,type="cairo")
        
      } else {
        
      }
      
      ### end of locuscompare ###
      }
    } # end of tissue loop
    
  } # end of gene loop
  
  # saving the results
  
  write.csv2(results_tot,paste0("/mnt/sda/gobemi01/ANGPTL3/coloc_res/",name_gwas,"_GWAS_top_hits_",tissue,".csv"),row.names=F)
  
