#%%
import pandas as pd
import altair as alt
import numpy as np
import hdbscan
import networkx as nx
import statsmodels.api as sm
import bioframe as bf
from multiprocessing import Pool
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
                          .assign(LFC=lambda df_:np.log10(df_.rate/(DNA_df.shape[0]/DNA_df.end.max()))))

hdb_cluster_summary_df = (hdb_cluster_summary_df
 .assign(fdr=lambda df_:sm.stats.multipletests(df_.pval,method='fdr_bh')[1]))
# %%
cl_depth_df = (pd.DataFrame.from_dict(nx.shortest_path_length(g, source=DNA_df.shape[0]),orient='index')
 .reset_index()
 .query('index > @DNA_df.shape[0] -1')
 .reset_index(drop=True)
 .rename(columns={
    'index':'cl',
    0:'level'
 }))
# %%
hdb_cluster_summary_df = hdb_cluster_summary_df.merge(cl_depth_df,how='inner')

#%%
hdb_cluster_summary_df = (hdb_cluster_summary_df.merge((egde_list
 .rename(columns={'child':"cl"}).loc[:,['cl','lambda_val','child_size']]),how='left'))
#%%%
read_parent = egde_list.query('child < @DNA_df.shape[0]').parent
noise_seed = hdb_cluster_summary_df.query("cl in  @read_parent").query('LFC < 0').cl
egde_list.query('child < @DNA_df.shape[0]').query('parent in @noise_seed')
# %%
peak_seed = hdb_cluster_summary_df.query("cl in  @read_parent").query('LFC > 0').cl.to_numpy()

# %%
def get_enriched_ancestry(i):
    tmp_ancestry = np.array(list(nx.ancestors(g,i)))
    tmp_ancestry = np.hstack((tmp_ancestry,i))
    tmp_df = hdb_cluster_summary_df.query('cl in @tmp_ancestry').loc[:,['cl','LFC','fdr','width','level']]

    tmp_df = tmp_df.sort_values('level',ascending=False).reset_index(drop=True)
    if tmp_df.LFC[0] < 0 :
        return tmp_df.iloc[[0],:].query('LFC > 0').cl.to_numpy()
    elif tmp_df.LFC[1] < 0 :
        return tmp_df.iloc[[0],:].query('LFC > 0').cl.to_numpy()

    elif sum(tmp_df.LFC > 0) < tmp_df.shape[0] :
        return tmp_df.iloc[list(range(0,tmp_df.LFC.where(tmp_df.LFC < 0).first_valid_index())),:].query('LFC > 0').cl.to_numpy()
    else:
        return tmp_df.query('LFC > 0').cl.to_numpy()
# %%
get_enriched_ancestry(6942)
# %%
num_processes = 4
#%%
# Execute the map function in parallel
with Pool(num_processes) as pool:
    results = pool.map(get_enriched_ancestry, range(0,DNA_df.shape[0]-1))

# %%
read_paths = np.concatenate(results).ravel()
# %%
unique_cl, cl_counts = np.unique(read_paths, return_counts=True)

# %%
ancestor_count = pd.DataFrame({
   'cl': unique_cl,
   'count': cl_counts
}).merge(hdb_cluster_summary_df.loc[:,['cl','level']],how='left')
# %%
count_summary = ancestor_count.groupby('level').agg(c=('count','sum')).reset_index()
(alt.Chart(count_summary
          )
 .mark_bar().encode(
    alt.X("level:Q", 
          bin=True),
    y=alt.Y(
        'c'
        ),
))
# %%
