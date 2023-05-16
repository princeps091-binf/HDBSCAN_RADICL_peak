#%%
import pandas as pd
import altair as alt
import numpy as np
import hdbscan
import networkx as nx
import bioframe as bf
import multiprocessing
import random
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
lambda_val = np.unique(egde_list.lambda_val.sort_values(ascending=True))[::-1][1:-1]

def cluster_nr(i):
    tmp_cl = (clusterer.single_linkage_tree_
                        .get_clusters(1/i,min_cluster_size=2))
    return (bf.cluster(DNA_df
                .assign(cl=tmp_cl)
                .query('cl > -1')
                .groupby('cl')
                .agg(start=('start','min'),end=('end','max'))
                .reset_index(drop=True)
                .assign(chrom=chromo)).shape[0])

# %%
with multiprocessing.Pool(processes=4) as pool:
        
        # Using map_async method to perform square operation on all numbers parallely
        result_cln = pool.map(cluster_nr ,lambda_val)        
lambda_cl_df = pd.DataFrame({'distance':(1/np.array(lambda_val)),
                             'cl':result_cln})
#%%
ch_line = (alt.Chart((lambda_cl_df
            .assign(lw=lambda df_:np.log10(df_.distance),
                    lr=lambda df_:np.log10(df_.cl)))
            )
.mark_line(
)
.encode(
    alt.X("lw:Q"),
    alt.Y('cl:Q')))

ch_line
#%%
lambda_peak = 1/lambda_cl_df.sort_values('cl',ascending=False).distance.iloc[0]
# %%
cl_set = egde_list.query('lambda_val >= @lambda_peak').query('child > @DNA_df.shape[0]').child.drop_duplicates().to_numpy()
# %%
seed_ancestor = np.array([len(list(set(nx.ancestors(g,i)) & set(cl_set))) for i in cl_set])
seed_children = np.array([len(list(set(nx.descendants(g,i)) & set(cl_set))) for i in cl_set])

# %%
cl_tree_context_df = pd.DataFrame({
      'cl':cl_set,
      'ancestor_count':seed_ancestor,
      'children_count':seed_children
})
#%%
candidate_peak_cl_df = (cl_tree_context_df
 .query('ancestor_count < 1')
 .query('children_count > 0')
 .sort_values('children_count'))

# %%
def get_cl_reads_idx(i):
    tmp_descendants = np.array(list(nx.descendants(g,i)))
    return tmp_descendants[tmp_descendants < DNA_df.shape[0]]
#%%
(DNA_df.iloc[get_cl_reads_idx(203206),:]
 .groupby('chrom')
 .agg(start = ('start',min),
      end= ('end',max),
      read_count= ('start','count'))
  .reset_index())

def get_cl_details(i):
    return(DNA_df.iloc[get_cl_reads_idx(i),:]
    .groupby('chrom')
    .agg(start = ('start',min),
        end= ('end',max),
        read_count= ('start','count'))
    .reset_index())

# %%
peak_read_idx = np.concatenate([get_cl_reads_idx(i) for i in candidate_peak_cl_df.cl.to_list()])
#%%
cl_summarydf = pd.concat([get_cl_details(i) for i in candidate_peak_cl_df.cl.to_list()])
# %%
alt.Chart(cl_summarydf
          .assign(width=lambda df_:df_.end-df_.start)
          .assign(lw=lambda df_:np.log10(df_.width))
          ).transform_window(
    ecdf="cume_dist()",
    sort=[{"field": "width"}],
).mark_line(
    interpolate="step-after"
).encode(
    x="width:Q",
    y="ecdf:Q"
)
#%%
alt.Chart(cl_summarydf
          .assign(lw=lambda df_:np.log10(df_.read_count))
          ).transform_window(
    ecdf="cume_dist()",
    sort=[{"field": "read_count"}],
).mark_line(
    interpolate="step-after"
).encode(
    x="read_count:Q",
    y="ecdf:Q"
)

# %%
(alt.Chart(dedup_radicl_df.iloc[peak_read_idx,:])
.mark_point(
    filled=True,
    size=1
)
.encode(
    alt.X("start:Q"),
    alt.Y('DNA_start:Q')))

# %%
random_sample = random.sample(range(0, DNA_df.shape[0]), len(peak_read_idx))
(alt.Chart(dedup_radicl_df.iloc[random_sample,:])
.mark_point(
    filled=True,
    size=1
)
.encode(
    alt.X("start:Q"),
    alt.Y('DNA_start:Q')))

# %%
