#%%
import pandas as pd
import altair as alt
import numpy as np
import bioframe as bf
from sklearn.neighbors import NearestNeighbors
import statsmodels.api as sm
from scipy.stats import norm
import hdbscan
import subprocess
import multiprocessing
alt.data_transformers.disable_max_rows()

#%%
chromo = 'chr19'
file_format = "BED"
ext_param= "25"
tmp_folder = "/home/vipink/Documents/FANTOM6/HDBSCAN_RADICL_peak/data/processed/"
black_list_file = "/home/vipink/Documents/FANTOM6/data/annotation/hg38-blacklist.v2.bed"

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
# load black list annotation
black_list_df=pd.read_csv(black_list_file,sep="\t",header=None)
black_list_df.columns = ['chrom','start','end','label']

# filter out black list regions
clean_RNA_df = bf.subtract(radicl_df,black_list_df)

# %%
(clean_RNA_df
 .query("strand == '+'")
 .to_csv(f"{tmp_folder}RADICL_RNA_{chromo}.bed",
         sep="\t",header=False,index=False))
#%%
extend_radicl_df = pd.concat([(clean_RNA_df
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
subprocess.run(["macs3",
                "pileup",
                "-f",file_format,
                "-i",f"{tmp_folder}RADICL_RNA_{chromo}.bed",
                "--extsize",ext_param,
                "-o",f"{tmp_folder}RADICL_RNA_pileup_{chromo}.bdg"],shell=False)

# %%
RNA_pile_bdg=pd.read_csv(f"{tmp_folder}RADICL_RNA_pileup_{chromo}.bdg",sep="\t",header=None)
RNA_pile_bdg = (RNA_pile_bdg
                .rename(columns={
                    0:'chrom',
                    1:'start',
                    2:'end',
                    3:'read_count'
                }))
# %%
(RNA_pile_bdg
 .assign(width=lambda df_:df_.end-df_.start)).query('read_count == 0').sort_values('width')
# %%
alt.Chart(bf.closest(RNA_pile_bdg.query('read_count > 0')).query('distance > 0').sort_values('distance')
 .assign(lw=lambda df_:np.log10(df_.distance))
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
signal_gap = bf.closest(RNA_pile_bdg.query('read_count > 0')).query('distance > 0').distance.mean()
neigh_agg = bf.closest(RNA_pile_bdg.query('read_count > 0')).query('distance > 0').query('distance < @signal_gap')
#%%
RNA_tag_cluster_df = (bf.cluster(extend_radicl_df,on=['strand'],min_dist=signal_gap)
 .groupby('cluster')
 .agg(chrom=('chrom','first'),
      start=('cluster_start','min'),
      end=('cluster_end','max'),
      strand=('strand','first'),
      size=('cluster','count'))
 .reset_index()
 .assign(width=lambda df_:df_.end - df_.start)
 .sort_values('width'))

tag_cluster_mean_width = RNA_tag_cluster_df.query('size > 1').width.mean()
# %%
subprocess.run(["macs3",
                "bdgpeakcall",
                "-i",f"{tmp_folder}RADICL_RNA_pileup_{chromo}.bdg",
                "-c","1",
                "-l",np.array2string(np.int64(np.ceil(tag_cluster_mean_width))),
                "-g",np.array2string(np.int64(np.ceil(signal_gap))),
                "--no-trackline",
                "-o",f"{tmp_folder}RADICL_DNA_{chromo}_peaks.bed"],shell=False)
#%%
RNA_peak_df=(pd.read_csv(f"{tmp_folder}RADICL_DNA_{chromo}_peaks.bed",sep="\t",header=None)
             .iloc[:,[0,1,2]]
             .rename(columns={
                 0:'chrom',
                 1:'start',
                 2:'end'
             }))

# %%
chr_read_rate = clean_RNA_df.shape[0]/((clean_RNA_df.query("strand == '+'")).end.max()-(clean_RNA_df.query("strand == '+'")).start.min())
#%%
(bf.count_overlaps(RNA_peak_df,(clean_RNA_df
 .query("strand == '+'")))
 .assign(width=lambda df_:df_.end-df_.start)
 .assign(peak_rate = lambda df_:df_.loc[:,'count']/df_.width)
 .assign(FC=lambda df_:df_.peak_rate/chr_read_rate))
# %%
