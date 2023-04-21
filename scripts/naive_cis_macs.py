#%%
import subprocess
import os
import re
import bioframe as bf
import pandas as pd
#%%
RADICL_file="/home/vipink/Documents/FANTOM6/HDBSCAN_RADICL_peak/data/processed/chr16_filter_df.csv"
bg_RADICL_file="/home/vipink/Documents/FANTOM6/HDBSCAN_RADICL_peak/data/processed/chr16_bg.bed"
black_list_file = "/home/vipink/Documents/FANTOM6/data/annotation/hg38-blacklist.v2.bed"
#%%
file_format = "BED"
ext_param= "147"
chromo = "chr16"
tmp_folder = "/home/vipink/Documents/FANTOM6/HDBSCAN_RADICL_peak/data/processed/"
#%%
# load the read tables (background and total)
radicl_df = pd.read_csv(RADICL_file,delimiter="\t")
dedup_radicl_df = (radicl_df.loc[:,['chrom','start','end','strand','DNA_start','DNA_end','DNA_strand']]
                   .drop_duplicates())
DNA_df = (dedup_radicl_df.loc[:,['chrom','DNA_start','DNA_end']]
    .rename(columns={'DNA_start':'start',
                    'DNA_end':'end'})
    .sort_values('start'))
bg_radicl_df = pd.read_csv(bg_RADICL_file,delimiter="\t")

# load black list annotation
black_list_df=pd.read_csv(black_list_file,sep="\t",header=None)
black_list_df.columns = ['chrom','start','end','label']

# filter out black list regions
clean_DNA_df = bf.subtract(DNA_df,black_list_df)
clean_bg_df = bf.subtract(bg_radicl_df,black_list_df)
#%%
(clean_DNA_df
 .to_csv(f"{tmp_folder}RADICL_DNA_{chromo}.bed",
         sep="\t",header=True,index=False))
(bg_radicl_df
 .to_csv(f"{tmp_folder}{chromo}_DNA_singleton.bed",
         sep="\t",header=True,index=False))

#%%
# MACS processing
subprocess.run(["macs3",
                "pileup",
                "-f",file_format,
                "-i",f"{tmp_folder}RADICL_DNA_{chromo}.bed",
                "--extsize",ext_param,
                "-o",f"{tmp_folder}RADICL_DNA_smooth_{chromo}.bdg"],shell=False)
subprocess.run(["macs3",
                "pileup",
                "-f",file_format,
                "-i",f"{tmp_folder}{chromo}_DNA_singleton.bed",
                "--extsize",ext_param,
                "-o",f"{tmp_folder}RADICL_DNA_bg_{chromo}.bdg"],shell=False)
subprocess.run(["macs3",
                "bdgcmp",
                "-t",f"{tmp_folder}RADICL_DNA_smooth_{chromo}.bdg",
                "-c",f"{tmp_folder}RADICL_DNA_bg_{chromo}.bdg",
                "-m","qpois",
                "-p","1",
                "-o",f"{tmp_folder}RADICL_DNA_qval_{chromo}.bdg"],shell=False)
subprocess.run(["macs3",
                "bdgpeakcall",
                "-i",f"{tmp_folder}RADICL_DNA_qval_{chromo}.bdg",
                "-c","1",
                "-l","75",
                "-g",ext_param,
                "--no-trackline",
                "-o",f"{tmp_folder}RADICL_DNA_{chromo}_peaks.bed"],shell=False)

# %%
