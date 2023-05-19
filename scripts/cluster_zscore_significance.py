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
RADICL_file="/home/vipink/Documents/FANTOM6/HDBSCAN_RADICL_peak/data/processed/chr16_filter_df.csv"
#%%
radicl_df = pd.read_csv(RADICL_file,delimiter="\t")

# %%
# remove duplicates
dedup_radicl_df = (radicl_df.loc[:,['chrom','start','end','strand','DNA_start','DNA_end','DNA_strand']]
                   .drop_duplicates()
                   .rename(columns={'start':'RNA_start',
                                    'end':'RNA_end',
                                    'strand':'RNA_strand'})
                   .sort_values('DNA_start')
                   .reset_index(drop=True)

)
# table for GAM
read_df = (dedup_radicl_df
            .loc[:,['RNA_start','DNA_start']]
            )
#%%
clusterer = hdbscan.HDBSCAN(min_cluster_size=2,
                            cluster_selection_epsilon=25,
                            metric='euclidean')
clusterer.fit(dedup_radicl_df.loc[:,['DNA_start']])
egde_list = clusterer.condensed_tree_.to_pandas()
g = clusterer.condensed_tree_.to_networkx()

#%%
nbrs = NearestNeighbors(n_neighbors=2, metric='euclidean',radius=25).fit(read_df.to_numpy())
distances, indices = nbrs.kneighbors(read_df.to_numpy())
# %%
read_neighbour_df = (dedup_radicl_df
 .assign(closest_DNA=np.abs(read_df.loc[indices[:,0],'DNA_start'].to_numpy() - read_df.loc[indices[:,1],'DNA_start'].to_numpy()),
         closest_RNA=np.abs(read_df.loc[indices[:,0],'RNA_start'].to_numpy() - read_df.loc[indices[:,1],'RNA_start'].to_numpy()))
         )
# %%
read_neighbour_df = (read_neighbour_df
 .assign(d=lambda df_: np.abs(df_.RNA_start - df_.DNA_start))
 )

# %%
read_neighbour_df = (read_neighbour_df
                     .assign(mod_RNA=lambda df_: df_.closest_RNA.where(df_.closest_RNA.gt(0),1),
                             mod_d = lambda df_: df_.d.where(df_.d.gt(0),1))
                     .assign(log_neigh=lambda df_:np.log10(df_.mod_RNA),
                             log_d=lambda df_:np.log10(df_.mod_d))
)
# %%
x_spline = read_neighbour_df[['log_d']]
y = read_neighbour_df.log_neigh.to_numpy()

bs = sm.gam.BSplines(x_spline, df=5, degree=3)

# %%
chr_gam = sm.GLMGam(y,smoother=bs)
chr_gam_res = chr_gam.fit()
# %%
gam_infl = chr_gam_res.get_influence()
# %%
read_neighbour_df = (read_neighbour_df
           .assign(zscore=gam_infl.resid_studentized)
            )
# %%
read_summary_df = (dedup_radicl_df
 .merge(read_neighbour_df,how='inner')
 .loc[:,['chrom','RNA_start','RNA_end','RNA_strand','DNA_start','DNA_end','DNA_strand','d','zscore']])
# %%
alt.Chart(read_summary_df).transform_density(
    'zscore',
    as_=['zscore', 'density'],
).mark_area().encode(
    x="zscore:Q",
    y='density:Q',
)

# %%
def get_cluster_zscore_signficance(i):
    tmp_descendant = np.array(list(nx.descendants(g,i)))
    tmp_read = tmp_descendant[tmp_descendant < read_summary_df.shape[0]]
    zsum = read_summary_df.iloc[tmp_read,:].zscore.sum()
    return norm(loc=0,scale=np.sqrt(len(tmp_read))).cdf(zsum)

# %%
hd_clusters = np.array(list(g.nodes))[np.array(list(g.nodes))>(read_summary_df.shape[0]-1)]

with multiprocessing.Pool(processes=4) as pool:
        # Using map_async method to perform square operation on all numbers parallely
        result_cln = pool.map(get_cluster_zscore_signficance ,hd_clusters)        

# %%
cl_zscore_significance_df = pd.DataFrame({
      'cl':hd_clusters,
      'pval':result_cln
}).sort_values('pval')
# %%
alt.Chart(cl_zscore_significance_df).mark_bar().encode(
    alt.X("pval:Q", bin=True),
    y='count()',
)
# %%
alt.Chart(cl_zscore_significance_df).transform_window(
    ecdf="cume_dist()",
    sort=[{"field": "pval"}],
).mark_line(
    interpolate="step-after"
).encode(
    x="pval:Q",
    y="ecdf:Q"
)

# %%
lambda_val = np.unique(egde_list.lambda_val.sort_values(ascending=True))[::-1][1:-1]

def cluster_nr(i):
    tmp_cl = (clusterer.single_linkage_tree_
                        .get_clusters(1/i,min_cluster_size=2))
    return sum(np.unique(tmp_cl) >= 0)
 
# %%
with multiprocessing.Pool(processes=5) as pool:
        # Using map_async method to perform square operation on all numbers parallely
        result_cln = pool.map(cluster_nr ,lambda_val)        
lambda_cl_df = pd.DataFrame({'distance':(1/np.array(lambda_val)),
                             'cl':result_cln})
# %%
lambda_cl_df
# %%
lambda_peak = 1/lambda_cl_df.sort_values('cl',ascending=False).distance.iloc[0]

# %%
cl_set = egde_list.query('lambda_val >= @lambda_peak').query('child > @dedup_radicl_df.shape[0]').child.drop_duplicates().to_numpy()

# %%
cl_zscore_significance_df.query('cl in @cl_set')
# %%

# %%
seed_ancestor = np.array([len(list(set(nx.ancestors(g,i)) & set(cl_set))) for i in cl_set])
seed_children = np.array([len(list(set(nx.descendants(g,i)) & set(cl_set))) for i in cl_set])
cl_tree_context_df = pd.DataFrame({
      'cl':cl_set,
      'ancestor_count':seed_ancestor,
      'children_count':seed_children
})

# %%
cl_zscore_significance_df = cl_zscore_significance_df.merge(cl_tree_context_df,how='inner')
#%%
cl_zscore_significance_df.query('ancestor_count == 0')
# %%
