#%%
import pandas as pd
import altair as alt
import numpy as np
import hdbscan
import networkx as nx
import statsmodels.api as sm
from joblib import Parallel, delayed

#%%
RADICL_read_file = "/home/vipink/Documents/FANTOM6/HDBSCAN_RADICL_peak/data/processed/chr16_filter_df.csv"

#%%
radicl_df = pd.read_csv(RADICL_read_file,delimiter="\t")
dedup_radicl_df = (radicl_df.loc[:,['chrom','start','end','strand','DNA_start','DNA_end','DNA_strand']]
                   .drop_duplicates())
DNA_df = (dedup_radicl_df.loc[:,['chrom','DNA_start','DNA_end']]
    .rename(columns={'DNA_start':'start',
                    'DNA_end':'end'})
    .sort_values('start'))

# %%
clusterer = hdbscan.HDBSCAN(min_cluster_size=2,
                            cluster_selection_epsilon=25,
                            metric='euclidean')
clusterer.fit(DNA_df.loc[:,['start']])

# %%
read_cluster_df = (DNA_df
           .assign(proba=clusterer.probabilities_,
                   outlier_score=clusterer.outlier_scores_,
                   out=clusterer.labels_ < 0,
                   cl=clusterer.labels_)
           .reset_index(drop=True))

#%%
g = clusterer.condensed_tree_.to_networkx()
egde_list = clusterer.condensed_tree_.to_pandas()

# %%
# get all cluster
hd_clusters = np.array(list(g.nodes))[np.array(list(g.nodes))>(DNA_df.shape[0]-1)]
#%%
tmp_cl = hd_clusters[58799]
# get descendants for specific cluster
tmp_cl_read_idx = np.array(list(nx.descendants(g,tmp_cl)))[np.array(list(nx.descendants(g,tmp_cl))) < DNA_df.shape[0]]
# get considered cluster reads
cl_summary_df = (read_cluster_df.iloc[tmp_cl_read_idx,:]
 .groupby('chrom')
 .agg(start = ('start','min'),
      end = ('end','max'),
      read_count = ('cl','count'))
 .reset_index()
 .assign(width=lambda df_:df_.end - df_.start)
 .assign(rate= lambda df_:df_.read_count/df_.width))
sm.stats.test_poisson_2indep(count1 = cl_summary_df.read_count[0],
                              exposure1 = cl_summary_df.width[0],
                              count2 = DNA_df.shape[0],
                              exposure2 = DNA_df.end.max(),
                              alternative='larger').pvalue
# %%
hdb_cluster_df = pd.DataFrame({
    'HDB_cluster': hd_clusters
})
# %%

def produce_cl_summary(name):
    tmp_cl = name
    # get descendants for specific cluster
    tmp_cl_read_idx = np.array(list(nx.descendants(g,tmp_cl)))[np.array(list(nx.descendants(g,tmp_cl))) < DNA_df.shape[0]]
    # get considered cluster reads
    return (read_cluster_df.iloc[tmp_cl_read_idx,:]
    .groupby('chrom')
    .agg(start = ('start','min'),
        end = ('end','max'),
        read_count = ('cl','count'))
    .reset_index()
    .assign(width=lambda df_:df_.end - df_.start)
    .assign(rate= lambda df_:df_.read_count/df_.width))

def produce_cl_read_idx(name):
    # get descendants for specific cluster
    tmp_cl_read_idx = np.array(list(nx.descendants(g,name)))[np.array(list(nx.descendants(g,name))) < DNA_df.shape[0]]
    return tmp_cl_read_idx

def produce_cl_start(idx_array):
    # get descendants for specific cluster
    return read_cluster_df.loc[idx_array,'start'].min()

def produce_cl_end(idx_array):
    # get descendants for specific cluster
    return read_cluster_df.loc[idx_array,'end'].max()

# %%
read_idx_list = list(map(produce_cl_read_idx, hdb_cluster_df.HDB_cluster.to_list()))

# %%
read_idx_list = Parallel(n_jobs=4)(delayed(produce_cl_read_idx)(i) for i in hdb_cluster_df.HDB_cluster.to_list())
# %%
tmp_cl_start = np.array(list(map(produce_cl_start, read_idx_list)))
tmp_cl_end = np.array(list(map(produce_cl_end, read_idx_list)))
tmp_cl_count = np.array(list(map(len, read_idx_list)))
# %%

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
# %%
hdb_cluster_summary_df = (hdb_cluster_summary_df
                          .assign(LFC=lambda df_:df_.rate/(DNA_df.shape[0]/DNA_df.end.max())))
