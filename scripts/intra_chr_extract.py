#%%
import ibis
import pandas as pd
import altair as alt

#%%
RADICL_read_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/Neuron/replicate1/raw/15.GGCTAC.2.bed"
black_list_file = "/home/vipink/Documents/FANTOM6/data/annotation/hg38-blacklist.v2.bed"
annotation_file = "/home/vipink/Documents/FANTOM6/data/annotation/FANTOM_CAT.lv3_robust.bed"

#%%
t = ibis.read_csv(RADICL_read_file,delim='\t')
# %%
t = t.relabel(dict(
    column00='RNA_chrom',
    column01='RNA_start',
    column02='RNA_end',
    column03='RNA_ID',
    column04='RNA_score',
    column05='RNA_strand',
    column06='DNA_chrom',
    column07='DNA_start',
    column08='DNA_end',
    column09='DNA_ID',
    column10='DNA_score',
    column11='DNA_strand',
    column12='RNA_label',
    column13='DNA_label',
    column14='inter_label'))
t.head().execute()
#%%
chromo = "chr19"
t_intra = t.filter((t.RNA_chrom == chromo) & (t.DNA_chrom == chromo)).execute()

# %%
alt.data_transformers.disable_max_rows()

(alt.Chart(t_intra)
.mark_point(
    size=0.1,
    filled=True,
    opacity=0.5
)
.encode(
    alt.X("RNA_start:Q"),
    alt.Y('DNA_start:Q')))

# %%
