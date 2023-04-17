#%%
import pandas as pd
import altair as alt
import bioframe as bf
import numpy as np
import hdbscan

#%%
RADICL_read_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/Neuron/replicate1/raw/RADICL_intra.bed"
black_list_file = "/home/vipink/Documents/FANTOM6/data/annotation/hg38-blacklist.v2.bed"
annotation_file = "/home/vipink/Documents/FANTOM6/data/annotation/FANTOM_CAT.lv3_robust.bed"

#%%
chromo = "chr16"

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
nascent_check_df = (transcript_read_inter
 .groupby(["chrom_1","RNA_ID_2","start_2","end_2","DNA_start_2","DNA_end_2"])
 .agg(edge_start=('start_1','min'),
      edge_end=('end_1','max'))
 .reset_index()
 .assign(max_d=lambda df_:np.abs(df_.start_2 - df_.edge_end),
         max_u=lambda df_:np.abs(df_.end_2 - df_.edge_start))
 .assign(max_e=lambda df_:df_[['max_d','max_u']].max(axis=1))
 .assign(max_r_d=lambda df_:np.abs(df_.start_2 - df_.DNA_end_2),
         max_r_u=lambda df_:np.abs(df_.end_2 - df_.DNA_start_2))
 .assign(max_r_e=lambda df_:df_[['max_r_d','max_r_u']].max(axis=1))
 .assign(nascent=lambda df_:df_.max_r_e < df_.max_e)
)
# %%
alt.data_transformers.disable_max_rows()

(alt.Chart(df
           .query("RNA_ID in @nascent_check_df.query('~nascent').RNA_ID_2"))
#           .query("strand == '+'"))
.mark_point(
    size=0.1,
    filled=True,
    opacity=0.5
)
.encode(
    alt.X("start:Q"),
    alt.Y('DNA_start:Q')))
#%%
(df
.query("RNA_ID in @nascent_check_df.query('~nascent').RNA_ID_2")
.to_csv(f"./../data/processed/{chromo}_filter_df.csv",
        sep="\t",
        header=True,index=False)
)

# %%

# %%
