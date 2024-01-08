#!/usr/bin/env Rscript

################################### #
#                                   #
#                                   #
#        ANGPTL3 project            #
#                                   #
#                                   #
################################### #

#### Information on exposure / outcome ####
# Exposure = TG UKB (no ANGPTL3 region, + top ANGPTL3 liver eQTL, + top ANGPTL3 blood pQTL)
# Outcomes = ApoB UKB 

top_pqtl = "rs10889352"
top_eqtl = "rs11207978"  

## Basic parameters ####

project = "ANGPTL3"

wd_data = "/mnt/sda/gobemi01/ANGPTL3"

wd_res_multi = "/mnt/sda/gobemi01/ANGPTL3/TG_ApoB"

setwd(wd_data)


## Load library ####

library(tidyverse)
library(data.table)
library(stringr)
library(TwoSampleMR)
library(coloc)
library(dplyr)
library(MendelianRandomization)
library(writexl)
library(LDlinkR)
library(ieugwasr)
library(Gobeil.R)

source("/mnt/sda/gobemi01/Functions/custom_functions.R")

window=1000000
threshold_pval_clumping = 5e-8
threshold_pval=5e-8
r2_clump=0.1

# ANGPTL3 region to remove #
mystart = 63063191
myend =  63071984

# Exposure ####

tg_gwas <- extract_instruments(outcomes="ieu-b-111",r2=0.1, p1=5e-08,clump=FALSE)
pqtl_tg <- tg_gwas %>% filter(SNP == top_pqtl)
eqtl_tg <- tg_gwas %>% filter(SNP == top_eqtl)
exp_dat <- extract_instruments(outcomes="ieu-b-111",r2=0.001, p1=5e-08,clump=TRUE)
to_remove <- exp_dat %>% filter(chr.exposure == "1" & between(pos.exposure, mystart-window/2, myend+window/2) == T)
exposure <- exp_dat %>% filter(SNP %in% to_remove$SNP == F )
exposure$exposure <- "TG_UKB"
exposure <- exposure %>% mutate(beta.exposure=beta.exposure*0.822)
exposure <- exposure %>% mutate(se.exposure=se.exposure*0.822)

# Outcome ####
pqtl_apob <- extract_outcome_data(snps=pqtl_tg$SNP,outcomes="ieu-b-108")
eqtl_apob <- extract_outcome_data(snps=eqtl_tg$SNP,outcomes="ieu-b-108")
outcome<-extract_outcome_data(snps=exposure$SNP,outcomes="ieu-b-108")
outcome$outcome <- "ApoB_UKB"
outcome <- outcome %>% mutate(beta.outcome=beta.outcome*0.24)
outcome <- outcome %>% mutate(se.outcome=se.outcome*0.24)
  # Harmonise data ####
  
  dat = harmonise_data(exposure, outcome, action=2)
  dat$exposure = "TG_UKB"
  dat$outcome = "ApoB_UKB"

dat_pqtl <- harmonise_data(pqtl_tg, pqtl_apob, action=2)
dat_pqtl$exposure = "TG_UKB_pQTL"
dat_pqtl$outcome = "ApoB_UKB_pQTL"
dat_eqtl <- harmonise_data(eqtl_tg, eqtl_apob, action=2)
dat_eqtl$exposure = "TG_UKB_eQTL"
dat_eqtl$outcome = "ApoB_UKB_eQTL"
  
res <- mr(dat,method_list = "mr_ivw")
dat <- dat %>% filter(mr_keep==TRUE)

res_pqtl <- mr(dat_pqtl,method_list = "mr_wald_ratio")
res_eqtl <- mr(dat_eqtl,method_list = "mr_wald_ratio") 


    #Label and create plots.
    library(ggrepel)
    library(ggplot2)
    source("/mnt/sda/gobemi01/Functions/custom_functions.R")
    
    ## SNP Label
    SNP_ANGPTL3<-c("rs10889352","rs11207978")
    dat_pqtl$label="rs10889352 (pQTL)"
    dat_eqtl$label="rs11207978 (eQTL)"
    dat <- dat %>% 
      mutate(label=ifelse(SNP %in% SNP_ANGPTL3, SNP,"")
      )
    index <- dat$beta.exposure > 0
    dat$beta.exposure[index] <- dat$beta.exposure[index] *-1
    dat$beta.outcome[index] <- dat$beta.outcome[index] *-1
    dat_tot <- rbind(dat,dat_pqtl,dat_eqtl)
  
    res_tot <- rbind(res,res_pqtl,res_eqtl)
    
    p1<-mr_scatter_plot_custom_custom(res,dat_tot,"ApoB_UKB",label_col = dat_tot$label)
    p1[[1]]
    
    #Save plot as a PDF.
    setwd(paste0("/mnt/sde/gobemi01/",project,"/TG_ApoB/"))
    ggsave(filename = paste0("/mnt/sde/gobemi01/",project,"/TG_ApoB/",project,"_convers_mr_TG_UKB_ApoB_UKB.png"),plot=last_plot(),width=7, height=7)
    
    fwrite(dat_tot,"/mnt/sda/gobemi01/ANGPTL3/TG_ApoB/data_convers_MR_Top_ANGPTL3_eQTL.rs11207978_pQTL.rs10889352_effect_TG_UKB_on_ApoB_UKB.csv",sep=";",dec=".")
    fwrite(res_tot,"/mnt/sda/gobemi01/ANGPTL3/TG_ApoB/res_convers_MR_Top_ANGPTL3_eQTL.rs11207978_pQTL.rs10889352_effect_TG_UKB_on_ApoB_UKB.csv",sep=";",dec=".")
    