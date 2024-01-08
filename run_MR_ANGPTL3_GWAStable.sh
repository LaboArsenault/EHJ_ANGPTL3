#!/usr/bin/env bash


# ************ 'Example of command to run the PHEWAS script' ************* #
# rs11207978 (top ANGPTL3 eQTL) # rs10889352 (top ANGPTL3 pQTL)
######################## RUN for GWAS_TABLE EXAMPLE ###########################################

nohup Rscript MR_ANGPTL3_singleSNP_GWAStable.R --wd /mnt/sda/gobemi01/PheWAS/PheWAS_MR_ANGPTL3/ --window=1000000 --pval=5e-08 --r2=0.1 \
--study gwas_table --gwas_table /mnt/sda/gobemi01/Data_instruments/gwas_table_13-03-2023.xlsx --genes_ref_file ./data/GRCh38_chr_pos_genes.txt \
--plink_binaries ./data/EUR_rs --r2_corr=0.6 --cis_trans cis --exp_start_pos 62597520 --exp_end_pos 62606313 \
--singleSNP_multiSNP singleSNP --split_var exposure --exp_pheno ANGPTL3_deCODE --exp_gwas_file=/mnt/sda/dbs_Web/deCODE/10391_1_ANGPTL3_ANGL3.txt.gz --exp_beta_col Beta \
--exp_se_col SE --exp_pval_col Pval --exp_chr_col Chrom --exp_eaf_col ImpMAF --exp_ea_col effectAllele --exp_oa_col otherAllele --exp_ss_col N \
--exp_snp_col rsids --exp_pos_col Pos --exp_p_is_log FALSE --exp_is_vcf=FALSE --exp_snp_col_num=4 --exp_chr_col_num=1 --exp_pos_col_num=2 \
--exp_a1_col_num=5 --exp_a2_col_num=6 --exp_maf_col_num=12 --exp_beta_col_num=7 --exp_se_col_num=10 --exp_p_col_num=8 --exp_n_col_num=11 \
--exp_ss_total=35559 --min_ratio=0.001 --min_cases=500 --to_exclude_file ./data/to_exclude_gwas_table.txt \
--out_p_is_log=FALSE --out_is_vcf=FALSE --save_every=1 --save_dat=TRUE --nthreads=10 \
--proxies=TRUE --reverse=FALSE --check_memory=TRUE --check_memory_lim=0.9 >> ./nohup/nohup.GWAStable_top_singleSNP_deCODE_w1000_r2_0.1.txt &