#%%
import pandas as pd
import altair as alt
import numpy as np
import bioframe as bf
from sklearn.neighbors import NearestNeighbors
import statsmodels.api as sm
from scipy.stats import norm
import hdbscan
import networkx as nx
import multiprocessing
alt.data_transformers.disable_max_rows()

#%%
chromo = 'chr1'
RNA_radicl = f'/home/vipink/Documents/FANTOM6/data/RADICL_data/IPSC/replicate1/raw/RNA/chr/RADICL_iPSC_RNA_{chromo}.bed'
#%%
radicl_df = pd.read_csv(RNA_radicl,delimiter="\t",header=None)
#%%
radicl_df = (radicl_df
             .rename(columns={
                 0:'chrom',
                 1:'start',
                 2:'end',
                 3:'ID',
                 4:'score',
                 5:'strand'
             }))
# %%
extend_radicl_df = pd.concat([(radicl_df
 .query("strand == '-'")
 .assign(start_b=lambda df_:df_.end - 25)
 .loc[:,['chrom','start_b','end','strand']]
 .rename(columns={'start_b':'start'})),
(radicl_df
 .query("strand == '+'")
 .assign(end_b=lambda df_:df_.start + 25)
 .loc[:,['chrom','start','end_b','strand']]
 .rename(columns={'end_b':'end'}))
]
)
plus_rate = extend_radicl_df.query("strand == '+'").shape[0]/(extend_radicl_df.query("strand == '+'").end.max() - extend_radicl_df.query("strand == '+'").start.min())
minus_rate = extend_radicl_df.query("strand == '-'").shape[0]/(extend_radicl_df.query("strand == '-'").end.max() - extend_radicl_df.query("strand == '-'").start.min())
rate_df = pd.DataFrame({
    'strand':['+','-'],
    'rate':[plus_rate,minus_rate]
})
#%%
RNA_tag_cluster_df = (bf.cluster(extend_radicl_df,on=['strand'])
 .groupby('cluster')
 .agg(chrom=('chrom','first'),
      start=('cluster_start','min'),
      end=('cluster_end','max'),
      strand=('strand','first'),
      size=('cluster','count'))
 .reset_index()
 .assign(width=lambda df_:df_.end - df_.start)
 .sort_values('width'))
#%%
plus_max_d = bf.closest(RNA_tag_cluster_df.query("strand == '+'").query('size > 1')).distance.mean()
minus_max_d = bf.closest(RNA_tag_cluster_df.query("strand == '-'").query('size > 1')).distance.mean()
max_d = max([plus_max_d,minus_max_d])

alt.Chart(bf.closest(RNA_tag_cluster_df.query("strand == '+'"))
          .assign(lw=lambda df_:np.log10(df_.distance))).transform_density(
    'lw',         # Specify the column containing the data
    as_=['agg_width', 'density']   # Name the new columns
).mark_area().encode(
    x='agg_width:Q', 
    y='density:Q')

#%%
RNA_clusters_df = (bf.cluster(extend_radicl_df,on=['strand'],min_dist=max_d)
 .groupby('cluster')
 .agg(chrom=('chrom','first'),
      start=('cluster_start','min'),
      end=('cluster_end','max'),
      strand=('strand','first'),
      size=('cluster','count'))
 .reset_index()
 .assign(width=lambda df_:df_.end - df_.start)
 .sort_values('width'))
#%%
candidate_enriched_cluster_df = RNA_clusters_df.merge(rate_df).assign(cluster_rate= lambda df_:df_.loc[:,'size']/df_.width).query('cluster_rate > rate')
#%%
candidate_enriched_cluster_df.loc[:,'size'].sum()/extend_radicl_df.shape[0]
#%%
candidate_enriched_cluster_df.query("strand == '+'").loc[:,'width'].sum()/(extend_radicl_df.query("strand == '+'").end.max() - extend_radicl_df.query("strand == '+'").start.min())

# %%
bar = alt.Chart(candidate_enriched_cluster_df.query("strand == '+'")
                .sort_values('start')

                ).mark_errorbar(
                    thickness=10
                ).encode(
    alt.X("end:Q",scale=alt.Scale(zero=False),title="coord(bp)"),
    alt.X2("start:Q"),
    alt.Y("size:Q"),
    strokeWidth=alt.value(1)
)

bar

# %%