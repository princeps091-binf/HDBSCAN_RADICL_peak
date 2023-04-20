#%%
import pandas as pd
import altair as alt
import numpy as np
from sklearn.neighbors import NearestNeighbors
import statsmodels.api as sm
from scipy.stats import t


alt.data_transformers.disable_max_rows()
#%%
RADICL_file="/home/vipink/Documents/FANTOM6/HDBSCAN_RADICL_peak/data/processed/chr16_filter_df.csv"
#%%
radicl_df = pd.read_csv(RADICL_file,delimiter="\t")

# %%
plus_strand_space_df = (radicl_df
                        .query("strand == '+'")
                        .loc[:,['start','DNA_start']]
                        .rename(columns={'start':'RNA_start'})
                        .reset_index(drop=True)
                        )
#%%
nbrs = NearestNeighbors(n_neighbors=2, metric='euclidean',radius=25).fit(plus_strand_space_df.to_numpy())
distances, indices = nbrs.kneighbors(plus_strand_space_df.to_numpy())
# %%
read_neighbour_df = (plus_strand_space_df
 .assign(closest_DNA=np.abs(plus_strand_space_df.loc[indices[:,0],'DNA_start'].to_numpy() - plus_strand_space_df.loc[indices[:,1],'DNA_start'].to_numpy()),
         closest_RNA=np.abs(plus_strand_space_df.loc[indices[:,0],'RNA_start'].to_numpy() - plus_strand_space_df.loc[indices[:,1],'RNA_start'].to_numpy())))
# %%
read_neighbour_df = (read_neighbour_df
 .assign(d=lambda df_: np.abs(df_.RNA_start - df_.DNA_start))
 .sort_values('d'))
# %%
chart = (alt.Chart(read_neighbour_df
           .assign(log_neigh=lambda df_:np.log10(df_.closest_DNA),
                   log_d=lambda df_:np.log10(df_.d))
)
.mark_point(opacity=0.5)
.encode(
    x="log_d:Q",
    y='log_neigh:Q'
))

chart + chart.transform_regression('log_d', 'log_neigh',method='linear').mark_line(color='red')

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
print(chr_gam_res.summary())
# %%
chr_gam_res.plot_partial(0, cpr=False)
# %%
gam_infl = chr_gam_res.get_influence()
# %%
gam_infl.resid_studentized
# %%
