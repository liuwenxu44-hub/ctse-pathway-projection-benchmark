"""No refits. Fixed paired summaries; contexts never pooled across scales."""
import csv,json,pathlib,sys
import numpy as np
import pandas as pd
def read_tsv(path):
 return pd.read_csv(path,sep='\t',keep_default_na=False,na_values=['NA'])
D=pathlib.Path(sys.argv[1]);out=D/'SUMMARIES';out.mkdir(exist_ok=False)
def save(df,n):df.to_csv(out/n,sep='\t',index=False,na_rep='NA',float_format='%.17g',lineterminator='\n')
unitkeys=['task_id','scenario','replicate_id','method','gene_stratum']
def balanced_unit(df,value,extra=()):
 keys=list(extra)+unitkeys
 bycell=df.groupby(keys+['celltype','group'],dropna=False,sort=True)[value].agg(['mean','count','size']).reset_index()
 # Both equally weighted groups are required when a cell/stratum has any
 # evaluable samples. Empty design strata are explicitly unavailable.
 a=bycell.groupby(keys+['celltype'],dropna=False,sort=True).agg(score=('mean',lambda x:x.mean() if x.notna().sum()==2 else np.nan),evaluable_samples=('count','sum'),planned_samples=('size','sum')).reset_index()
 u=a.groupby(keys,dropna=False,sort=True).agg(score=('score','mean'),evaluable_celltypes=('score','count'),planned_celltypes=('score','size'),evaluable_samples=('evaluable_samples','sum'),planned_samples=('planned_samples','sum')).reset_index()
 return a,u
def distribution(u,extra=()):
 keys=list(extra)+['scenario','method','gene_stratum']
 return u.groupby(keys,sort=True,dropna=False).agg(planned_units=('score','size'),evaluable_units=('score','count'),mean=('score','mean'),median=('score','median'),q25=('score',lambda x:x.quantile(.25)),q75=('score',lambda x:x.quantile(.75))).reset_index()
rank=read_tsv(D/'GENE_RANK_RECOVERY_LONG.tsv')
a,u=balanced_unit(rank,'spearman');save(a,'RANK_CELLTYPE_GROUPBALANCED.tsv');save(u,'RANK_UNIT_SUMMARY.tsv');save(distribution(u),'RANK_SCENARIO_DISTRIBUTION.tsv')
join=['task_id','scenario','replicate_id','celltype','group','sample_id','gene_stratum']
pairs=[]
for b in ['reference','bulk']:
 z=rank.merge(rank[rank.method==b][join+['spearman']],on=join,how='left',suffixes=('','_baseline'),validate='many_to_one')
 z=z[z.method!=b].copy();z['difference']=z.spearman-z.spearman_baseline;_,pu=balanced_unit(z,'difference');pu['baseline']=b;pu['difference_definition']='method_minus_baseline_Spearman';pairs.append(pu)
paired=pd.concat(pairs,ignore_index=True);save(paired,'RANK_PAIRED_BASELINE_UNIT.tsv');save(distribution(paired,extra=['baseline','difference_definition']),'RANK_PAIRED_BASELINE_SCENARIO.tsv')
del rank
err=read_tsv(D/'SAME_SCALE_ERROR_LONG.tsv');errunit=[]
for metric in ['rmse','sne']:
 _,eu=balanced_unit(err,metric,extra=['context']);eu['metric']=metric;errunit.append(eu)
eu=pd.concat(errunit,ignore_index=True);save(eu,'SAME_SCALE_ERROR_UNIT.tsv');save(distribution(eu,extra=['context','metric']),'SAME_SCALE_ERROR_SCENARIO.tsv')
del err
dr=read_tsv(D/'SCALE_SPECIFIC_DIRECTION_LONG.tsv');summary=[];conf=[];recalls=[]
classes=['negative','zero','positive']
def direction_group(context,method,path,scope,z):
 truth=z.truth_direction;pred=z.estimated_direction;ok=pred.notna();valid_truth=truth.isin(classes);ok=ok&valid_truth
 recall=[]
 for cl in classes:
  take=truth==cl;n=int((take&ok).sum());correct=int((take&ok&(pred==cl)).sum());val=correct/n if n else np.nan;recall.append(val)
  recalls.append(dict(context=context,method=method,pathway=path,scope=scope,truth_class=cl,planned=int(take.sum()),evaluable=n,correct=correct,recall=val,coverage=n/int(take.sum()) if take.sum() else np.nan))
  for pc in classes+['NO_OUTPUT']:
   num=int((take&((~ok) if pc=='NO_OUTPUT' else (ok&(pred==pc)))).sum());conf.append(dict(context=context,method=method,pathway=path,scope=scope,truth=cl,prediction=pc,n=num))
 nz=truth.isin(['negative','positive']);zero=truth=='zero'
 summary.append(dict(context=context,method=method,pathway=path,scope=scope,independent_unit='simulation_replicate',planned_replicates=z.replicate_id.nunique() if scope!='ALL_SCENARIOS' else z.task_id.nunique(),planned_records=len(z),valid_truth_records=int(valid_truth.sum()),evaluable_records=int(ok.sum()),coverage=ok.mean(),overall_accuracy=(truth[ok]==pred[ok]).mean() if ok.any() else np.nan,balanced_accuracy_macro_recall=np.mean(recall) if all(np.isfinite(recall)) else np.nan,nonzero_planned=int(nz.sum()),nonzero_evaluable=int((nz&ok).sum()),nonzero_sign_accuracy=(truth[nz&ok]==pred[nz&ok]).mean() if (nz&ok).any() else np.nan,zero_planned=int(zero.sum()),zero_evaluable=int((zero&ok).sum()),zero_truth_recall=(pred[zero&ok]=='zero').mean() if (zero&ok).any() else np.nan))
for (ctx,m,path),z in dr.groupby(['context','method','pathway'],sort=True):
 direction_group(ctx,m,path,'ALL_SCENARIOS',z)
 for sc,v in z.groupby('scenario',sort=True):direction_group(ctx,m,path,sc,v)
save(pd.DataFrame(summary),'DIRECTION_ACCURACY_AND_COVERAGE.tsv');save(pd.DataFrame(conf),'THREE_CLASS_CONFUSION_MATRIX.tsv');save(pd.DataFrame(recalls),'CLASS_RECALL_AND_COVERAGE.tsv')
contrast=read_tsv(D/'SCALE_SPECIFIC_CONTRAST_LONG.tsv');keys=['context','scenario','method','gene_stratum','celltype_design_changed']
# Cell types are first averaged within simulation unit, never treated as independent replicates.
cu=contrast.groupby(keys+['task_id','replicate_id'],sort=True,dropna=False).agg(contrast_error_L2=('contrast_error_L2','mean'),estimated_contrast_L2=('estimated_contrast_L2','mean'),truth_contrast_L2=('truth_contrast_L2','mean'),evaluable_celltypes=('contrast_error_L2','count')).reset_index();save(cu,'CONTRAST_UNIT.tsv')
cs=cu.groupby(keys,sort=True,dropna=False).agg(planned_units=('task_id','size'),evaluable_units=('contrast_error_L2','count'),mean_contrast_error_L2=('contrast_error_L2','mean'),median_contrast_error_L2=('contrast_error_L2','median'),median_estimated_contrast_L2=('estimated_contrast_L2','median'),median_truth_contrast_L2=('truth_contrast_L2','median')).reset_index();save(cs,'CONTRAST_SCENARIO.tsv')
availability=read_tsv(D/'METHOD_AVAILABILITY.tsv');av=availability.groupby(['scenario','method'],sort=True).agg(planned_units=('task_id','size'),available_units=('available','sum'),min_returned=('returned_genes','min'),max_returned=('returned_genes','max'),min_common=('common_genes','min'),max_common=('common_genes','max')).reset_index();av['coverage']=av.available_units/av.planned_units;save(av,'METHOD_AVAILABILITY_SUMMARY.tsv')
(out/'SUMMARY_COMPLETE').write_text('Deterministic, paired and scale-separated; no model refits.\n')
print('SUMMARY_COMPLETE',len(summary),'direction summary rows')
