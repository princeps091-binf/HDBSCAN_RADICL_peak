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
 .to_csv(f"{tmp_folder}RADICL_RNA_{chromo}.bed",
         sep="\t",header=False,index=False))

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
alt.Chart(RNA_pile_bdg
 .assign(width=lambda df_:df_.end-df_.start).query('read_count == 0')
 .assign(lw=lambda df_:np.log10(df_.width))
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
