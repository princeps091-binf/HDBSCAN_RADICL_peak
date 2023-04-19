#%%
import pandas as pd
import bioframe as bf
alt.data_transformers.disable_max_rows()
#%%
RADICL_file="/home/vipink/Documents/FANTOM6/HDBSCAN_RADICL_peak/data/processed/chr16_filter_df.csv"
#%%
radicl_df = pd.read_csv(RADICL_file,delimiter="\t")
# Observe likely PCR duplicates (tiny minority)
#%%
dedup_radicl_df = (radicl_df.loc[:,['chrom','start','end','strand','DNA_start','DNA_end','DNA_strand']]
                   .drop_duplicates())

# %%
def produce_bg_read_set(radicl_df):
    cluster_read_df = (bf.cluster(radicl_df))

    cluster_size_df = (cluster_read_df
    .groupby('cluster')
    .agg(size=('cluster','count'))
    .reset_index()
    .sort_values('size'))

    cluster_size_summary = (cluster_size_df
    .groupby('size')
    .agg(cl_s=('size','sum'))
    .reset_index()
    .assign(cums=lambda df_:df_.cl_s.cumsum()/radicl_df.shape[0]))


    cluster_size_thresh = cluster_size_summary.loc[:,'size'].loc[cluster_size_summary.cums.gt(0.5).idxmax()]

    bg_cluster_set = (cluster_size_df
    .query("size <= @cluster_size_thresh")
    .loc[:,'cluster'])

    return(cluster_read_df
    .query('cluster in @bg_cluster_set')
    .sort_values('start')
    .loc[:,['chrom','start','end']])
# %%
radicl_df = (dedup_radicl_df.loc[:,['chrom','DNA_start','DNA_end']]
    .rename(columns={'DNA_start':'start',
                    'DNA_end':'end'})
    .sort_values('start'))
bg_read_df = produce_bg_read_set(radicl_df)
#%%
bf.closest(bg_read_df).sort_values('distance').agg(med=('distance','median'))
# %%
