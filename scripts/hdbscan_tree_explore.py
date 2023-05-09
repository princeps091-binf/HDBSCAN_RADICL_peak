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

#%%
hdb_cluster_summary_df.to_csv("./../data/processed/chr16_hdb_cl_summary.csv",
                              sep="\t",
                              index=False,
                              header=True)
#%%
local_cluster_df = (hdb_cluster_summary_df
 .query('width < 10000')
 .query('LFC > 0')
 .sort_values('fdr'))

#%%
def make_cluster_thresh(thresh):
    return (bf.cluster(local_cluster_df
            .query('fdr <= @thresh'))
            .groupby(['cluster','cluster_start','cluster_end'])
            .agg(ncl=('cluster','count'))
            .reset_index()
            .assign(width=lambda df_:df_.cluster_end-df_.cluster_start,
                    fdr=thresh).shape[0])

def make_cluster_thresh_w(thresh):
    return (bf.cluster(local_cluster_df
            .query('fdr <= @thresh'))
            .groupby(['cluster','cluster_start','cluster_end'])
            .agg(ncl=('cluster','count'))
            .reset_index()
            .assign(width=lambda df_:df_.cluster_end-df_.cluster_start,
                    fdr=thresh)
            .agg(scl=('width','sum')).width.to_numpy())

pval_thresh = local_cluster_df.fdr.unique()
#%%
# Define the number of processes to use
num_processes = 4
#%%
# Execute the map function in parallel
with Pool(num_processes) as pool:
    results_w = pool.map(make_cluster_thresh_w, pval_thresh)
    results_n = pool.map(make_cluster_thresh, pval_thresh)
#%%
res_df = pd.DataFrame({
    'fdr':pval_thresh,
    'ncl':np.hstack(results_w),
    'wcl':np.hstack(results_n)
})
(alt.Chart(res_df
           .assign(lf=lambda df_:np.log10(df_.fdr),
                   lncl=lambda df_:np.log10(df_.ncl),
                   lw=lambda df_:np.log10(df_.wcl))
           .assign(clp=lambda df_:df_.lf.where(df_.lf >= -50,-50)))
.mark_line().encode(
    x="clp:Q",
    y="lncl:Q"
))
#%%
zero_bump = hdb_cluster_summary_df.pval.to_numpy()[hdb_cluster_summary_df.pval.to_numpy()> 0].min()/2

(alt.Chart(hdb_cluster_summary_df
            .query('width < 10000')
            .assign(lw=lambda df_:np.log10(df_.width),
                    lr=lambda df_:np.log10(df_.read_count),
                    lpval=lambda df_:-np.log10(df_.pval + zero_bump),
                    qp=lambda df_:pd.qcut(-np.log10(df_.fdr + zero_bump),10,labels=pd.qcut(-np.log10(df_.fdr + zero_bump),10,retbins=True)[1][0:-1])
                    )
            .assign(clp=lambda df_:df_.lpval.where(df_.lpval <= 20,20))
            )
.transform_window(
    ecdf="cume_dist()",
    sort=[{"field": "clp"}],
).mark_line(
    interpolate="step-after"
).encode(
    x="clp:Q",
    y="ecdf:Q"
))
#%%
(alt.Chart(hdb_cluster_summary_df
 .assign(lpval=lambda df_:-np.log10(df_.pval + zero_bump))
 .sort_values('pval')
 .assign(clp= lambda df_:df_.lpval.cumsum()/df_.lpval.sum(),
         perc_rank=np.array(list(range(0,hdb_cluster_summary_df.shape[0])))/hdb_cluster_summary_df.shape[0],
         rank = np.array(list(range(0,hdb_cluster_summary_df.shape[0]))))
 )
 .mark_line(
    interpolate="step-after"
).encode(
    x="rank:Q",
    y="clp:Q"
)
)
#%%
bar = alt.Chart(hdb_cluster_summary_df
                .sort_values('pval')
                .assign(
                        perc_rank=np.array(list(range(0,hdb_cluster_summary_df.shape[0])))
                )
                .iloc[0:900000,:]
                .query('LFC > 0')
                .query('start > 50200000')
                .query('end < 50500000')
                .query('width < 10000')
                .query('pval < 10**(-10)')
                ).mark_errorbar().encode(
    alt.X("end:Q",scale=alt.Scale(zero=False),title="coord(bp)"),
    alt.X2("start:Q"),
    alt.Y("level:O"),
    alt.Color("LFC:Q",scale=alt.Scale(scheme='plasma',reverse=False))
)

bar
#%%
(alt.Chart(hdb_cluster_summary_df
            .assign(lw=lambda df_:np.log10(df_.width),
                    lr=lambda df_:np.log10(df_.read_count),
                    lpval=lambda df_:-np.log10(df_.fdr + zero_bump))
            )
.mark_point(
     filled=True,
     opacity=0.2
)
.encode(
    alt.X("LFC:Q"),
    alt.Y('lpval:Q'),
    alt.Color('lw',scale = alt.Scale(scheme='viridis'))
    ))

#%%
zero_bump = hdb_cluster_summary_df.pval.to_numpy()[hdb_cluster_summary_df.pval.to_numpy()> 0].min()/2

(alt.Chart(hdb_cluster_summary_df
            .assign(lw=lambda df_:np.log10(df_.width),
                    lr=lambda df_:np.log10(1/df_.lambda_val),
                    lpval=lambda df_:-np.log10(df_.fdr + zero_bump))
            .assign(qp=lambda df_:pd.qcut(-np.log10(df_.fdr + zero_bump),10,labels=pd.qcut(-np.log10(df_.fdr + zero_bump),10,retbins=True)[1][0:-1]))
            )
.mark_point(
     filled=True,
     opacity=1
)
.encode(
    alt.X("lw:Q"),
    alt.Y('LFC:Q'),
    alt.Color('qp:O',scale = alt.Scale(scheme='blueorange'))
    ))

# %%
read_ancestry_df = pd.DataFrame({'ancestry':([np.sort(np.array(list(nx.ancestors(g,i)))) for i in range(0,DNA_df.shape[0])])})
#%%
(read_ancestry_df
 .assign(ancestor_n=lambda df_:df_.ancestry.map(len))
 .sort_values('ancestor_n'))


# %%

(alt.Chart(read_ancestry_df)
.mark_bar()
.encode(
    alt.X('ancestry:Q'),
    alt.Y('count()'),
))
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
