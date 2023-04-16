#%%
import pandas as pd
import altair as alt
import bioframe as bf
from joblib import Parallel, delayed
import multiprocessing
import numpy as np

#%%
RADICL_read_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/Neuron/replicate1/raw/RADICL_intra.bed"
black_list_file = "/home/vipink/Documents/FANTOM6/data/annotation/hg38-blacklist.v2.bed"
annotation_file = "/home/vipink/Documents/FANTOM6/data/annotation/FANTOM_CAT.lv3_robust.bed"

#%%
chromo = "chr22"

iter_csv = pd.read_csv(RADICL_read_file, iterator=True, chunksize=100000,sep='\t',header=None)

df = pd.concat([chunk.rename(columns={0:'chrom'}).query('chrom == @chromo') for chunk in iter_csv])
# %%
df = df.rename(columns={
    'chrom':'RNA_chrom',
    1:'RNA_start',
    2:'RNA_end',
    3:'RNA_ID',
    4:'RNA_score',
    5:'RNA_strand',
    6:'DNA_chrom',
    7:'DNA_start',
    8:'DNA_end',
    9:'DNA_ID',
    10:'DNA_score',
    11:'DNA_strand',
    12:'RNA_label',
    13:'DNA_label',
    14:'inter_label'})
#%%
transcript_annotation_df = pd.read_csv(annotation_file,header=None,delimiter="\t")

transcript_annotation_df.columns = ['chrom','start','end',
                                    'ID','score','strand',
                                    'start_b','end_b','sym',
                                    'exon_count','exon_length','exon_start']

#%%
df = df.rename(columns={'RNA_chrom':'chrom',
                   'RNA_start':'start',
                   'RNA_end':'end',
                   'RNA_strand':'strand'})
transcript_read_inter = bf.overlap(transcript_annotation_df.loc[:,['chrom','start','end','strand','ID']], df, 
                                       on=['strand'],
                                       how='inner', 
                                       suffixes=('_1','_2'))

#%%
def check_nascent(df_,name):
    return(df_
     .assign(start_d=lambda df_: np.abs(df_.start_2 - df_.start_1),
         end_d=lambda df_: np.abs(df_.end_2 - df_.end_1))
     .assign(max_d= lambda df_: df_.loc[:,['start_d','end_d']].max(axis=1))
     .groupby('RNA_ID_2')
     .agg(top_d=('max_d',max),
          end=('DNA_end_2',lambda x: x.iloc[0]),
          start=('DNA_start_2',lambda x: x.iloc[0]),
          RNA_start=('start_2',lambda x: x.iloc[0]))
     .reset_index()
     .assign(start_d=lambda x:np.abs(x.RNA_start - x.start)<x.top_d,
             end_d =lambda x: np.abs(x.RNA_start - x.end)<x.top_d)
     .assign(nascent = lambda x: sum([x.start_d,x.end_d])>0)
     .reset_index()
     .loc[:,['RNA_ID_2','nascent']])

def applyParallel(dfGrouped, func):
    retLst = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(func)(group,name) for name, group in dfGrouped)
    return pd.concat(retLst)

# %%
nascent_check_df = applyParallel(transcript_read_inter .groupby("RNA_ID_2"), check_nascent)
#%%
nascent_check_df.query('~nascent').RNA_ID_2
# %%
alt.data_transformers.disable_max_rows()

(alt.Chart(df.query("RNA_ID in @nascent_check_df.query('nascent').RNA_ID_2"))
.mark_point(
    size=0.1,
    filled=True,
    opacity=0.5
)
.encode(
    alt.X("start:Q"),
    alt.Y('DNA_start:Q')))

# %%
