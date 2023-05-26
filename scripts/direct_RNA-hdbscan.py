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
plus_strand_agg_df = RNA_tag_cluster_df.query("size > 1 & strand == '+'").assign(mid = lambda df_:df_.start + np.ceil(df_.width/2))
# %%

clusterer = hdbscan.HDBSCAN(min_cluster_size=2,
                            cluster_selection_epsilon=25,
                            metric='euclidean')
clusterer.fit(plus_strand_agg_df.loc[:,['mid']])

#%%
(plus_strand_agg_df
 .assign(cluster= clusterer.labels_)
 .query('cluster < 0')
)

# %%
(plus_strand_agg_df
 .assign(cluster= clusterer.labels_)
 .query('cluster > -1')
 .groupby('cluster')
 .agg(agg_count=('start','count'))
 .reset_index()
 .sort_values('agg_count'))
# %%
bf.count_overlaps(plus_strand_agg_df
 .assign(cluster= clusterer.labels_)
 .query('cluster == 6559')
 .groupby('cluster')
 .agg(
     chrom=('chrom','first'),
     start=('start','min'),
     end=('end','max'),
     strand=('strand','first'))
 .reset_index(),extend_radicl_df,on=['strand'])

# %%
bf.count_overlaps((plus_strand_agg_df
 .assign(cluster= clusterer.labels_)
 .query('cluster < 0')
),extend_radicl_df,on=['strand']
)
# %%
