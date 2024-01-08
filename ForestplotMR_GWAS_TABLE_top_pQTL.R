####################### ANGPTL3 pQTL MR FORESTPLOT #############################

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(ggplot2)
library(MRInstruments)
library(forestplot)
library(Rmpfr)

wd_data = "/mnt/sda/gobemi01/ANGPTL3/PheWAS/Results_IEU_GWAS/MultiSNP_MR/"

project = "ANGPTL3"

### Import all MR results ###

setwd(wd_data)

gwas_res <- fread("/mnt/sda/gobemi01/PheWAS/PheWAS_MR_ANGPTL3/results/MRres/gwas_table/ANGPTL3_deCODE/w1000_r2_0.1/MR_results_ANGPTL3_deCODE_cis_singleSNP_pval.5e-08_w1000_r2.0.1.txt")

apob_res <- gwas_res[1,]
ldl_res <- gwas_res[14,]
hdl_res <- gwas_res[10,]
tg_res <- gwas_res[22,]
tg_res$abs_zscore = abs(tg_res$"b"/tg_res$"se")
a <- 2*Rmpfr::pnorm(q=mpfr(tg_res$abs_zscore,precBits = 53),lower.tail = FALSE,log.p = FALSE)
tg_res$pval <- formatMpfr(a,3)
cad_res <- gwas_res[8,]
hf_res <- gwas_res[11,]
is_res <- gwas_res[12,] 
ap_res <- gwas_res[3,]   
ckd_res <- gwas_res[4,]
nafld_res <- gwas_res[17,]
t2d_res <- gwas_res[23,]   
  
##### Format betas ####

data_forest <- matrix(nrow = 15, ncol=3)
colnames(data_forest) = c("mean", "lower", "upper")
rownames(data_forest) <- c('Outcomes','Lipid and lipoprotein levels' ,'Apo B','HDL','LDL','TG',
                           'Cardiovascular diseases','Coronary artery disease','Heart failure','Ischemic stroke',
                           'Metabolic diseases','Acute pancreatitis','Chronic kidney disease','Non-alcoholic fatty liver disease','Type 2 diabetes')
data_forest[1,] <- NA
data_forest[2,] <- NA
data_forest[3,] <- c(apob_res$b, apob_res$b_CIlow,apob_res$b_CIhigh)
data_forest[4,] <- c(hdl_res$b, hdl_res$b_CIlow,hdl_res$b_CIhigh)
data_forest[5,] <- c(ldl_res$b, ldl_res$b_CIlow,ldl_res$b_CIhigh)
data_forest[6,] <- c(tg_res$b, tg_res$b_CIlow,tg_res$b_CIhigh)
data_forest[7,] <- NA
data_forest[8,] <- c(cad_res$b, cad_res$b_CIlow,cad_res$b_CIhigh)
data_forest[9,] <- c(hf_res$b, hf_res$b_CIlow,hf_res$b_CIhigh)
data_forest[10,] <- c(is_res$b, is_res$b_CIlow,is_res$b_CIhigh)
data_forest[11,] <- NA
data_forest[12,] <- c(ap_res$b, ap_res$b_CIlow,ap_res$b_CIhigh)
data_forest[13,] <- c(ckd_res$b, ckd_res$b_CIlow,ckd_res$b_CIhigh)
data_forest[14,] <- c(nafld_res$b, nafld_res$b_CIlow,nafld_res$b_CIhigh)
data_forest[15,] <- c(t2d_res$b, t2d_res$b_CIlow,t2d_res$b_CIhigh)

## LABELS

label <- matrix(nrow=15, ncol=4)

label[1,1] <- 'Outcomes'
label[2,1] <- 'Lipid and lipoprotein levels'
label[3,1] <- 'Apo B'
label[4,1] <- 'HDL-C'
label[5,1] <- 'LDL-C'
label[6,1] <- 'TG'
label[7,1] <- 'Cardiovascular diseases'
label[8,1] <- 'Coronary artery disease'
label[9,1] <- 'Heart failure'
label[10,1] <- 'Ischemic stroke'
label[11,1] <- 'Metabolic diseases'
label[12,1] <- 'Acute pancreatitis'
label[13,1] <- 'Chronic kidney disease'
label[14,1] <- 'Non-alcoholic fatty liver disease'
label[15,1] <- 'Type 2 diabetes'

format_OR_SE = function(x) {format(round(x,digits = 2), nsmall = 2)}
label[1,3] <- 'Beta (95% CI)'
label[2,3] <- NA
label[3,3] <- paste0(format_OR_SE(apob_res$b),' (', format_OR_SE(apob_res$b_CIlow),'-',format_OR_SE(apob_res$b_CIhigh),')')
label[4,3] <- paste0(format_OR_SE(hdl_res$b),' (', format_OR_SE(hdl_res$b_CIlow),'-',format_OR_SE(hdl_res$b_CIhigh),')')
label[5,3] <- paste0(format_OR_SE(ldl_res$b),' (', format_OR_SE(ldl_res$b_CIlow),'-',format_OR_SE(ldl_res$b_CIhigh),')')
label[6,3] <- paste0(format_OR_SE(tg_res$b),' (', format_OR_SE(tg_res$b_CIlow),'-',format_OR_SE(tg_res$b_CIhigh),')')
label[7,3] <- NA
label[8,3] <- paste0(format_OR_SE(cad_res$b),' (', format_OR_SE(cad_res$b_CIlow),'-',format_OR_SE(cad_res$b_CIhigh),')')
label[9,3] <- paste0(format_OR_SE(hf_res$b),' (', format_OR_SE(hf_res$b_CIlow),'-',format_OR_SE(hf_res$b_CIhigh),')')
label[10,3] <- paste0(format_OR_SE(is_res$b),' (', format_OR_SE(is_res$b_CIlow),'-',format_OR_SE(is_res$b_CIhigh),')')
label[11,3] <- NA
label[12,3] <- paste0(format_OR_SE(ap_res$b),' (', format_OR_SE(ap_res$b_CIlow),'-',format_OR_SE(ap_res$b_CIhigh),')')
label[13,3] <- paste0(format_OR_SE(ckd_res$b),' (', format_OR_SE(ckd_res$b_CIlow),'-',format_OR_SE(ckd_res$b_CIhigh),')')
label[14,3] <- paste0(format_OR_SE(nafld_res$b),' (', format_OR_SE(nafld_res$b_CIlow),'-',format_OR_SE(nafld_res$b_CIhigh),')')
label[15,3] <- paste0(format_OR_SE(t2d_res$b),' (', format_OR_SE(t2d_res$b_CIlow),'-',format_OR_SE(t2d_res$b_CIhigh),')')

format_pval_sci = function(x) {format(x, nsmall = 3,scientific=TRUE,digits = 3)}
format_pval_dec = function(x) {format(x, nsmall = 3,scientific=FALSE,digits = 2)}
label[1,4] <- 'P-Value'
label[2,4] <- NA
label[3,4] <- format_pval_sci(apob_res$pval)
label[4,4] <- format_pval_sci(hdl_res$pval)
label[5,4] <- format_pval_sci(ldl_res$pval)
label[6,4] <- tg_res$pval
label[7,4] <- NA
label[8,4] <- format_pval_dec(cad_res$pval)
label[9,4] <- format_pval_dec(hf_res$pval)
label[10,4] <- format_pval_dec(is_res$pval)
label[11,4] <- NA
label[12,4] <- format_pval_dec(ap_res$pval)
label[13,4] <- format_pval_dec(ckd_res$pval)
label[14,4] <- format_pval_dec(nafld_res$pval)
label[15,4] <- format_pval_dec(t2d_res$pval)

parametres_text = fpTxtGp(label = gpar(extrafont = "Arial", cex = 1), xlab  = gpar(fontfamily = "Arial", cex = 1), ticks = gpar(fontfamily = "", cex = 1.0))

forestplot::forestplot(
  label,
  mean = data_forest[,"mean"],
  lower = data_forest[,"lower"],
  upper = data_forest[,"upper"],
  is.summary = c(T,T,F,F,F,F,T,F,F,F,T,F,F,F,F),
  zero = 0,
  graphwidth = unit(x = 7, units = "cm"),
  xticks = c(seq(from = -0.2, to = 0.4, by = 0.1)),
  clip = c(0),
  graph.pos = 3,
  lwd.ci = 4,
  lwd.xaxis = 4,
  lwd.zero = 4,
  boxsize = 0.22,
  vertices=T,
  txt_gp = parametres_text,
  title="",
  xlab= substitute(paste("Effect size (95% CI) of genetically-predicted plasma ANGPTL3 levels on each outcome")))

# Save forestplot

