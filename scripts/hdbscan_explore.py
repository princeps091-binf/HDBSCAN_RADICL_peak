#%%
import pandas as pd
import altair as alt
import bioframe as bf
import numpy as np
import hdbscan
import networkx as nx

#%%
RADICL_file="/home/vipink/Documents/FANTOM6/HDBSCAN_RADICL_peak/data/processed/chr16_filter_df.csv"
#%%
radicl_df = pd.read_csv(RADICL_file,delimiter="\t")

# %%
plus_strand_space_df = (radicl_df
                        .query("strand == '+'")
                        .loc[:,['start','DNA_start']]
                        )

# %%
clusterer = hdbscan.HDBSCAN(min_cluster_size=2,
                            cluster_selection_epsilon=25,
                            metric='chebyshev',
                            cluster_selection_method='leaf')
clusterer.fit(plus_strand_space_df)

#%%
g = clusterer.condensed_tree_.to_networkx()
egde_list = clusterer.condensed_tree_.to_pandas()
#%%
nx.ancestors(g,)
# %%
hdbscan_read_edges = (egde_list
 .query('child_size == 1')
 .query('parent > @plus_strand_space_df.shape[0]')
 .sort_values('child'))
hdbscan_read_edges.groupby('parent').agg(size=('child','count')).reset_index().sort_values('size')

# %%
read_cluster_df = (plus_strand_space_df
           .assign(proba=clusterer.probabilities_,
                   outlier_score=clusterer.outlier_scores_,
                   out=clusterer.labels_ < 0,
                   cl=clusterer.labels_)
           .reset_index(drop=True)
           .assign(lambda_val=hdbscan_read_edges.reset_index().lambda_val))
#%%
exemplars = clusterer.exemplars_

# %%
alt.data_transformers.disable_max_rows()

(alt.Chart(read_cluster_df.assign(log_val=lambda df_:np.log10(df_.outlier_score + 1)))
.transform_density(
    'log_val',
    groupby=['out'],
    as_=['log_val', 'density']
).mark_area(opacity=0.3).encode(
    x="log_val:Q",
    y='density:Q',
    color='out'
))
# %%
cluster_span_df = (read_cluster_df
 .groupby('cl')
 .agg(RNA_min=('start','min'),
      RNA_max=('start','max'),
      DNA_min=('DNA_start','min'),
      DNA_max=('DNA_start','max'),
      size=('start','count'))
 .reset_index()
 .assign(RNA_span=lambda df_:df_.RNA_max - df_.RNA_min,
         DNA_span=lambda df_:df_.DNA_max - df_.DNA_min)
 .query('cl > -1')
)
# %%
alt.data_transformers.disable_max_rows()

(alt.Chart(cluster_span_df.assign(log_val=lambda df_:np.log10(df_.DNA_span)))
.transform_density(
    'log_val',
    as_=['log_val', 'density']
).mark_area(opacity=0.3).encode(
    x="log_val:Q",
    y='density:Q'
))

# %%
(read_cluster_df
.assign(cl_col=lambda df_:df_.cl.where(df_.cl.eq(0),-1)))
# %%
alt.data_transformers.disable_max_rows()

(alt.Chart((read_cluster_df
.assign(cl_col=lambda df_:df_.cl.where(df_.cl.isin(cluster_span_df.query('size > 10').cl),-1))))

.mark_point(
    size=0.1,
    filled=True,
    opacity=1
)
.encode(
    alt.X("start:Q"),
    alt.Y('DNA_start:Q'),
    color=alt.Color('cl_col:N', legend=None)))

# %%
nx.shortest_path_length(g,plus_strand_space_df.shape[0])
alt.data_transformers.disable_max_rows()
read_idx = np.array(list(nx.descendants(g,92969)))[np.array(list(nx.descendants(g,92969))) < plus_strand_space_df.shape[0]]
(alt.Chart((read_cluster_df
.iloc[read_idx,:]))

.mark_point(
    size=1,
    filled=True,
    opacity=1
)
.encode(
    alt.X("start:Q"),
    alt.Y('DNA_start:Q')))

# %%
nx.shortest_path_length(g,plus_strand_space_df.shape[0])
# %%
