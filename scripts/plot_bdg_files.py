#%%
import pandas as pd
import altair as alt
import numpy as np
import subprocess
alt.data_transformers.disable_max_rows()

#%%
chromo = 'chr16'
RADICL_file="/home/vipink/Documents/FANTOM6/HDBSCAN_RADICL_peak/data/processed/chr16_filter_df.csv"
file_format = "BED"
smooth= "1000"
tmp_folder = "/home/vipink/Documents/FANTOM6/HDBSCAN_RADICL_peak/data/processed/"

#%%
subprocess.run(["macs3",
                "pileup",
                "-f",file_format,
                "-i",f"{tmp_folder}RADICL_DNA_{chromo}.bed",
                "--extsize",smooth,
                "-o",f"{tmp_folder}RADICL_DNA_smooth_{smooth}_{chromo}.bdg"],shell=False)

bdg_file = f"{tmp_folder}RADICL_DNA_smooth_{smooth}_{chromo}.bdg"
#%%
bdg_df = (pd.read_csv(bdg_file,delimiter="\t",header=None)
          .drop(0)
          .rename(columns={
              0:'chrom',
              1:'start',
              2:'end',
              3:"value"
          }))

#%%
bar = alt.Chart(bdg_df
                .sort_values('start')
                .query('start>20000000')
                .query('end < 20005000')
                ).mark_errorbar().encode(
    alt.X("end:Q",scale=alt.Scale(zero=False),title="coord(bp)"),
    alt.X2("start:Q"),
    alt.Y("value:Q"),
    strokeWidth=alt.value(1)
)

bar

# %%
