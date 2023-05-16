#%%
import pandas as pd
import altair as alt
import numpy as np
import hdbscan
import networkx as nx
import statsmodels.api as sm
import bioframe as bf
import multiprocessing
alt.data_transformers.disable_max_rows()
#%%
chromo = 'chr16'
RADICL_read_file = f"/home/vipink/Documents/FANTOM6/HDBSCAN_RADICL_peak/data/processed/{chromo}_filter_df.csv"

#%%
radicl_df = pd.read_csv(RADICL_read_file,delimiter="\t")
dedup_radicl_df = (radicl_df.loc[:,['chrom','start','end','strand','DNA_start','DNA_end','DNA_strand']]
                   .drop_duplicates())
DNA_df = (dedup_radicl_df.loc[:,['chrom','DNA_start','DNA_end']]
    .rename(columns={'DNA_start':'start',
                    'DNA_end':'end'})
    .sort_values('start')
    .reset_index(drop=True))

# %%
clusterer = hdbscan.HDBSCAN(min_cluster_size=2,
                            cluster_selection_epsilon=25,
                            metric='euclidean')
clusterer.fit(DNA_df.loc[:,['start']])
#%%
egde_list = clusterer.condensed_tree_.to_pandas()
g = clusterer.condensed_tree_.to_networkx()
#%%
alt.Chart(egde_list
          .assign(lw=lambda df_:np.log10(1/df_.lambda_val))).transform_window(
    ecdf="cume_dist()",
    sort=[{"field": "lw"}],
).mark_line(
    interpolate="step-after"
).encode(
    x="lw:Q",
    y="ecdf:Q"
)
#%%
def compute_noise_level(i):
     return sum(clusterer.single_linkage_tree_.get_clusters(1/i, min_cluster_size=2) < 0)/tot_read_count

def lambda_cluster_size(i):
    return np.median(np.unique(clusterer.single_linkage_tree_.get_clusters(1/i, min_cluster_size=2),return_counts=True)[1][1:-1])

lambda_val = np.unique(egde_list.lambda_val.sort_values(ascending=True))[::-1]
tot_read_count = DNA_df.shape[0]
#%%
with multiprocessing.Pool(processes=4) as pool:
        result = pool.map(lambda_cluster_size ,lambda_val)        
#%%
lambda_noise_df = pd.DataFrame({'distance':(1/np.array(lambda_val)),'prop':result})

(alt.Chart((lambda_noise_df
            .assign(lw=lambda df_:np.log10(df_.distance),
                    lr=lambda df_:np.log10(df_.prop)))
            )
.mark_line(
)
.encode(
    alt.X("lw:Q"),
    alt.Y('lr:Q')))


# %%
# %%
