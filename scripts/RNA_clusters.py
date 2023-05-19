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
chromo = 'chr19'
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
RNA_cluster_df = (bf.cluster(extend_radicl_df,on=['strand'])
 .groupby('cluster')
 .agg(size=('cluster','count'),
      start=('cluster_start','min'),
      end=('cluster_end','max'),
      strand=('strand','first'))
 .reset_index()
 .assign(width=lambda df_:df_.end - df_.start)
 .sort_values('width'))
# %%
alt.Chart(RNA_cluster_df
          .assign(lw=lambda df_:np.log10(df_.width))
          ).transform_window(
    ecdf="cume_dist()",
    sort=[{"field": "lw"}],
).mark_line(
    interpolate="step-after"
).encode(
    x="lw:Q",
    y="ecdf:Q"
)

# %%
alt.Chart(RNA_cluster_df
          .sort_values('size',ascending=False)
          .assign(lw=lambda df_:np.log10(df_.width),
                  cs = lambda df_:df_.loc[:,'size'].cumsum()
                  )
          ).mark_line(
    interpolate="step-after"
).encode(
    x="size:Q",
    y="cs:Q"
)

# %%
bar = alt.Chart(RNA_cluster_df.query('size > 1').query("strand == '+'")
                .sort_values('start')
                .query('start>10000000')
                .query('end < 10200000')

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