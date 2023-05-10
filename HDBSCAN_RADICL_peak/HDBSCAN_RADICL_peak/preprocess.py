def detect_nascent_reads(transcript_annotation_df,df):
    transcript_read_inter = bf.overlap(transcript_annotation_df.loc[:,['chrom','start','end','strand','ID']], df, 
                                       on=['strand'],
                                       how='inner', 
                                       suffixes=('_1','_2'))


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
    return nascent_check_df
