#%%
import pandas as pd
import altair as alt
import numpy as np
from sklearn.neighbors import NearestNeighbors
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
plus_strand_space_df = (dedup_radicl_df
                        .query("strand == '+'")
                        .loc[:,['start','DNA_start']]
                        .reset_index()
                        )
#%%
nbrs = NearestNeighbors(n_neighbors=2, metric='euclidean',radius=25).fit(plus_strand_space_df.to_numpy())
distances, indices = nbrs.kneighbors(plus_strand_space_df.to_numpy())
# %%
read_neighbour_df = (plus_strand_space_df
 .assign(closest_DNA=np.abs(plus_strand_space_df.loc[indices[:,0],'DNA_start'].to_numpy() - plus_strand_space_df.loc[indices[:,1],'DNA_start'].to_numpy()),
         closest_RNA=np.abs(plus_strand_space_df.loc[indices[:,0],'start'].to_numpy() - plus_strand_space_df.loc[indices[:,1],'start'].to_numpy())))
# %%
dna_dist_cdf = (read_neighbour_df
 .sort_values('closest_DNA')
 .groupby('closest_DNA')
 .agg(read_count=('start','count'))
 .reset_index()
 .assign(cread=lambda df_:df_.read_count.cumsum()/plus_strand_space_df.shape[0])
 .rename(columns={'closest_DNA':'distance'})
 .assign(end='DNA'))

rna_dist_cdf = (read_neighbour_df
 .sort_values('closest_RNA')
 .groupby('closest_RNA')
 .agg(read_count=('start','count'))
 .reset_index()
 .assign(cread=lambda df_:df_.read_count.cumsum()/plus_strand_space_df.shape[0])
 .rename(columns={'closest_RNA':'distance'})
 .assign(end='RNA'))

tot_df = pd.concat([rna_dist_cdf,dna_dist_cdf])

# %%
(alt.Chart(tot_df.assign(log_val=lambda df_:np.log10(df_.distance + 1)))
.mark_line(opacity=0.6).encode(
    x="log_val:Q",
    y='cread:Q',
    color="end"
))
# %%
