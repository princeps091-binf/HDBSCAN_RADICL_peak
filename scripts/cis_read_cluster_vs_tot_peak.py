#%%
import bioframe as bf
import pandas as pd
import numpy as np
#%%
RADICL_read_file = "/home/vipink/Documents/FANTOM6/HDBSCAN_RADICL_peak/data/processed/chr16_filter_df.csv"
MACS_bdg_file = "/home/vipink/Documents/FANTOM6/HDBSCAN_RADICL_peak/data/processed/RADICL_DNA_qval_chr16.bdg"
tot_read_peak_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/IPSC/replicate1/results_clean/MACS/peaks/RADICL_DNA_chr16_peaks.bed"
#%%
radicl_df = pd.read_csv(RADICL_read_file,delimiter="\t")
dedup_radicl_df = (radicl_df.loc[:,['chrom','start','end','strand','DNA_start','DNA_end','DNA_strand']]
                   .drop_duplicates())
DNA_df = (dedup_radicl_df.loc[:,['chrom','DNA_start','DNA_end']]
    .rename(columns={'DNA_start':'start',
                    'DNA_end':'end'})
    .sort_values('start'))

MACS_bdg_df = pd.read_csv(MACS_bdg_file,delimiter="\t",header=None)
# %%
DNA_cluster_df = (bf.cluster(DNA_df)
                 .groupby(['cluster','chrom','cluster_start','cluster_end'])
                 .agg(nread=('start','count'))
                 .reset_index()
                 .rename(columns={'cluster_start':"start",
                                  "cluster_end":"end"}))
MACS_bdg_df = (MACS_bdg_df
               .rename(columns={0:"chrom",
                                1:"start",
                                2:"end",
                                3:"qval"}))
# %%
thresh = -np.log10(0.1)

candidate_peak = (bf.overlap(DNA_cluster_df,MACS_bdg_df,return_overlap=True)
 .query("~(chrom_ in [None])")
 .assign(qval_s=lambda df_:(df_.overlap_end - df_.overlap_start)*df_.qval_)
 .groupby(['cluster','nread','chrom','start','end'])
 .agg(qval=('qval_','mean'),
      nover=('start_','count'),
      qval_m=('qval_s','sum'))
 .reset_index()
 .assign(qval_cl=lambda df_:df_.qval_m/(df_.end-df_.start),
         cl_w=lambda df_:df_.end-df_.start)
 .query('qval > @thresh')
 .sort_values('qval_cl'))
# %%
tot_peak_df = pd.read_csv(tot_read_peak_file,delimiter="\t",header=None)
tot_peak_df = (tot_peak_df
               .rename(columns={
                    0:"chrom",
                    1:"start",
                    2:"end"}))
# %%
bf.closest(candidate_peak,tot_peak_df.loc[:,['chrom','start','end']])
# %%
