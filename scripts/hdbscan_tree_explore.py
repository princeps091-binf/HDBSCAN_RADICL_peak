#%%
import pandas as pd
import altair as alt
import numpy as np
import hdbscan
import networkx as nx
import statsmodels.api as sm
from joblib import Parallel, delayed
import bioframe as bf
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
g = clusterer.condensed_tree_.to_networkx()
egde_list = clusterer.condensed_tree_.to_pandas()
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
#Set distance thresh
clusterer.single_linkage_tree_.get_clusters(10000, min_cluster_size=2)
np.median(np.unique(clusterer.single_linkage_tree_.get_clusters(10000, min_cluster_size=2)[(clusterer.single_linkage_tree_.get_clusters(10000, min_cluster_size=2) > 0)],return_counts=True)[1])



# %%
# get all cluster
hd_clusters = np.array(list(g.nodes))[np.array(list(g.nodes))>(DNA_df.shape[0]-1)]
hdb_cluster_df = pd.DataFrame({
    'HDB_cluster': hd_clusters
})
#%%
def produce_cl_read_idx(name):
    # get descendants for specific cluster
    tmp_cl_read_idx = np.array(list(nx.descendants(g,name)))[np.array(list(nx.descendants(g,name))) < DNA_df.shape[0]]
    return tmp_cl_read_idx

def produce_cl_start(idx_array):
    # get descendants for specific cluster
    return DNA_df.loc[idx_array,'start'].min()

def produce_cl_end(idx_array):
    # get descendants for specific cluster
    return DNA_df.loc[idx_array,'end'].max()

# %%
read_idx_list = list(map(produce_cl_read_idx, hdb_cluster_df.HDB_cluster.to_list()))
tmp_cl_start = np.array(list(map(produce_cl_start, read_idx_list)))
tmp_cl_end = np.array(list(map(produce_cl_end, read_idx_list)))
tmp_cl_count = np.array(list(map(len, read_idx_list)))
#%%
hdb_cluster_summary_df=(pd.DataFrame({
    'cl': hdb_cluster_df.HDB_cluster,
    'start':tmp_cl_start,
    'end':tmp_cl_end,
    'read_count':tmp_cl_count})
    .assign(chrom='chr16')
    .assign(width=lambda df_:df_.end-df_.start)
    .assign(rate=lambda df_:df_.read_count/df_.width))
# %%
def check_rate_sign(df_,i):
    return (sm.stats.test_poisson_2indep(count1 = df_.read_count[i],
                                exposure1 = df_.width[i],
                                count2 = DNA_df.shape[0],
                                exposure2 = DNA_df.end.max(),
                                alternative='larger').pvalue)

# %%
tmp_pval = list(map(lambda i: check_rate_sign(hdb_cluster_summary_df,i), range(0,hdb_cluster_summary_df.shape[0])))
# %%
hdb_cluster_summary_df = hdb_cluster_summary_df.assign(pval=tmp_pval)
hdb_cluster_summary_df = (hdb_cluster_summary_df
                          .assign(LFC=lambda df_:np.log10(df_.rate/(1/np.mean(bf.closest(DNA_df).distance)))))

egde_list.assign(out=lambda df_:np.isfinite(df_.lambda_val)).query('~out')

#%%
alt.Chart(hdb_cluster_summary_df
          .assign(lw=lambda df_:np.log10(df_.width))).transform_window(
    ecdf="cume_dist()",
    sort=[{"field": "lw"}],
).mark_line(
    interpolate="step-after"
).encode(
    x="lw:Q",
    y="ecdf:Q"
)
#%%

(alt.Chart(hdb_cluster_summary_df
          .assign(lp=lambda df_:-np.log10(df_.pval))
          )
.transform_window(
    ecdf="cume_dist()",
    sort=[{"field": "level"}],
).mark_line(
    interpolate="step-after"
).encode(
    x="level:Q",
    y="ecdf:Q"
))
#%%
size_thresh = 1000
tmp_hdb_cluster_summary_df = (hdb_cluster_summary_df
                              .query('width <= @size_thresh'))

[nx.ancestors(g,i) for i in tmp_hdb_cluster_summary_df.cl]
# %%
read_ancestry_df = pd.DataFrame({'ancestry':np.array([len(nx.ancestors(g,i)) for i in range(0,DNA_df.shape[0])])})
# %%

(alt.Chart(read_ancestry_df)
.mark_bar()
.encode(
    alt.X('ancestry:Q'),
    alt.Y('count()'),
))
# %%
read_ancestry_df = pd.DataFrame({'ancestry':np.array([nx.ancestors(g,i) for i in range(0,DNA_df.shape[0])])})
#%%
def get_ancestry_stat(df_,g,i):
    read_ancestor = nx.ancestors(g,i)
    return df_.query('cl in @read_ancestor').pval.min()
#%%
ancestry_size_df = pd.DataFrame({'ancestry':np.hstack([get_ancestry_stat(hdb_cluster_summary_df,g,i) for i in range(0,123)])})

# %%
cl_depth_df = (pd.DataFrame.from_dict(nx.shortest_path_length(g, source=158032),orient='index')
 .reset_index()
 .query('index > @DNA_df.shape[0] -1')
 .reset_index(drop=True)
 .rename(columns={
    'index':'cl',
    0:'level'
 }))
# %%
hdb_cluster_summary_df = hdb_cluster_summary_df.merge(cl_depth_df,how='inner')
# %%
(alt.Chart((hdb_cluster_summary_df
            .assign(lp=lambda df_:-np.log10(df_.pval + 1),
                    lw=lambda df_:np.log10(df_.width))
            .assign(score=lambda df_:((df_.lp-df_.lp.mean())/df_.lp.std())* ((df_.LFC-df_.LFC.mean())/df_.LFC.std()))
            )
            )
.mark_point(
    filled=True,
    size=12,
    opacity=1
)
.encode(
    alt.X("lw:Q"),
    alt.Y('score:Q')))
# %%
(alt.Chart(hdb_cluster_summary_df
          .query('LFC > 0')
          .assign(lw=lambda df_:np.log10(df_.width),
                  lp=lambda df_:-np.log10(df_.pval))
          .query('lp<50'))
.transform_density(
    'lw',         # Specify the column containing the data
    as_=['agg_width', 'density']   # Name the new columns
).mark_area(
    opacity=0.3
).encode(
    x='agg_width:Q', 
    y='density:Q'))

# %%
sm.stats.test_p
# %%
