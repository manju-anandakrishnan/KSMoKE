import pandas as pd
import os

substrate_max_score_dfs = []
for feather_file in os.listdir('ksmo_predictions'):
    df = pd.read_feather('ksmo_predictions/'+feather_file)
    substrate_pred_max = df.groupby(['substrate_motif'])['ksf2_pred'].max()
    substrate_pred_max = substrate_pred_max.reset_index()
    substrate_max_score_dfs.append(substrate_pred_max)
    

substrate_max_score_df = pd.concat(substrate_max_score_dfs)
print(substrate_max_score_df.shape)

sm_df = pd.read_csv('substrates_motif_v1.csv',sep='|')
sm_df['sm'] = sm_df['Seq_Substrate']+'_'+sm_df['Motif']
merged_df = sm_df.merge(substrate_max_score_df,how='inner',left_on=['sm'],right_on=['substrate_motif'])
merged_df.drop(['sm','substrate_motif'],axis=1,inplace=True)
merged_df.rename({'ksf2_pred':'max_score'},axis=1,inplace=True)
merged_df.to_csv('substrate_motif.csv',sep='|',index=False)

