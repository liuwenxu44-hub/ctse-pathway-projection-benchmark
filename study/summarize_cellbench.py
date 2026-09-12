"""Deterministic endpoint aggregation; no model execution or outcome selection."""
import csv,sys,math,statistics as st,hashlib
from pathlib import Path
from collections import defaultdict,Counter
P=Path(sys.argv[1])
def read(name):
 with open(P/name,encoding='utf-8') as f:return list(csv.DictReader(f,delimiter='\t'))
def num(x):
 try:return float(x)
 except (TypeError,ValueError):return math.nan
def valid(x):return math.isfinite(num(x))
def mean(xs):
 xs=[num(x) for x in xs if valid(x)];return st.mean(xs) if xs else math.nan
def quantile(xs,p):
 xs=sorted(num(x) for x in xs if valid(x))
 if not xs:return math.nan
 i=(len(xs)-1)*p;j=int(i);return xs[j]+(xs[min(j+1,len(xs)-1)]-xs[j])*(i-j)
def write(name,rows):
 assert rows,name
 with open(P/name,'w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader()
  for r in rows:w.writerow({k:('NA' if isinstance(v,float) and not math.isfinite(v) else format(v,'.17g') if isinstance(v,float) else v) for k,v in r.items()})
def group(rows,keys):
 d=defaultdict(list)
 for r in rows:d[tuple(r[k] for k in keys)].append(r)
 return d
def aggregate(rows,keys,prefix):
 srows=[];crows=[];summary=[]
 for key,rr in sorted(group(rows,keys).items()):
  base=dict(zip(keys,key));vals=[];eligible_comps=0
  for co,cc in sorted(group(rr,['composition']).items()):
   ss=[];eligible=sum(str(r.get('eligible','True')).lower()=='true' for r in cc)
   for si in ['mRNA_3.75','mRNA_7.50','mRNA_15.00','mRNA_30.00']:
    child=[r for r in cc if r['stratum']==si];v=mean(r['value'] for r in child)
    n=sum(valid(r['value']) for r in child)
    srows.append(dict(**base,composition=co[0],stratum=si,n_planned=len(child),n_eligible=sum(str(r.get('eligible','True')).lower()=='true' for r in child),n_available=n,value=v,status='AVAILABLE' if valid(v) else 'STRATUM_ENDPOINT_UNAVAILABLE'))
    ss.append(v)
   v=mean(ss);vals.append(v);eligible_comps+=eligible>0
   crows.append(dict(**base,composition=co[0],n_strata_planned=4,n_strata_available=sum(valid(z) for z in ss),n_libraries_planned=len(cc),n_libraries_eligible=eligible,n_libraries_available=sum(valid(r['value']) for r in cc),value=v,status='AVAILABLE' if valid(v) else 'MIXTURE_ENDPOINT_UNAVAILABLE'))
  summary.append(dict(**base,n_compositions_planned=len(vals),n_compositions_eligible=eligible_comps,n_compositions_available=sum(valid(z) for z in vals),median=quantile(vals,.5),Q1=quantile(vals,.25),Q3=quantile(vals,.75),IQR=quantile(vals,.75)-quantile(vals,.25),cell_line_bias_mean=mean(vals) if base.get('metric')=='BIAS_EST_MINUS_ORACLE' else math.nan,independent_unit='DESIGN_COMPOSITION_NOT_DONOR'))
 write(prefix+'_STRATUM.tsv',srows);write(prefix+'_COMPOSITION.tsv',crows);write(prefix+'_SUMMARY.tsv',summary)
ct=read('CTSE_LIBRARY_LONG.tsv');fr=read('FRACTION_LIBRARY_LONG.tsv')
keys=['task','method','platform','arm','context','proxy','cell_line','metric']
assert len(group(ct,keys+['composition','stratum','library']))==len(ct)
aggregate(ct,keys,'CTSE')
# Exact inherited fraction aggregation; the scale-correction release had reused
# these unchanged tables. This public rebuild materializes them from that same
# frozen library-level table, without another fraction estimate.
aggregate(fr,['platform','method','arm','cell_line','metric'],'FRACTION')
# Paired comparisons are formed before scalar aggregation, never from mismatched native scores.
basenames=['REFERENCE_ONLY_CTSE','BULK_COPY_CTSE']
bk=['platform','context','proxy','cell_line','composition','stratum','library','metric']
indexes={m:{tuple(r[k] for k in bk):r for r in ct if r['method']==m and not r['context'].endswith('_NATIVE')} for m in basenames}
paired=[]
def pair(a,b,comparison,kind):
 assert a['gene_count']==b['gene_count']
 ok=valid(a['value']) and valid(b['value'])
 return dict(comparison=comparison,comparison_type=kind,method=a['method'],platform=a['platform'],arm=a['arm'],context=a['context'],proxy=a['proxy'],cell_line=a['cell_line'],composition=a['composition'],stratum=a['stratum'],library=a['library'],metric=a['metric'],gene_count=a['gene_count'],eligible=str(a['eligible']).lower()=='true' and str(b['eligible']).lower()=='true',left_value=num(a['value']),right_value=num(b['value']),left_status=a['status'],right_status=b['status'],value=num(a['value'])-num(b['value']) if ok else math.nan,status='AVAILABLE' if ok else 'PAIRED_ENDPOINT_UNAVAILABLE')
for a in ct:
 if a['method'] in basenames or a['context'].endswith('_NATIVE'):continue
 for baseline in basenames:
  b=indexes[baseline].get(tuple(a[k] for k in bk));assert b is not None
  paired.append(pair(a,b,a['method']+'_MINUS_'+baseline,'METHOD_MINUS_BASELINE'))
oracle={tuple(r[k] for k in ['method']+bk):r for r in ct if r['arm']=='ORACLE_FRACTION_SENSITIVITY' and '_PAIR_' in r['context']}
for a in ct:
 if a['arm']!='PRACTICAL_NONORACLE' or '_PAIR_' not in a['context']:continue
 b=oracle.get(tuple(a[k] for k in ['method']+bk));assert b is not None
 paired.append(pair(a,b,a['method']+'_PRACTICAL_MINUS_ORACLE','PRACTICAL_MINUS_ORACLE'))
pk=['task','method','platform','arm','context','cell_line','composition','stratum','library','metric']
c2={tuple(r[k] for k in pk):r for r in ct if r['proxy']=='C2'}
for a in ct:
 if a['proxy']=='C1':paired.append(pair(a,c2[tuple(a[k] for k in pk)],a['method']+'_C1_MINUS_C2','C1_MINUS_C2'))
write('PAIRED_LIBRARY_LONG.tsv',paired)
aggregate(paired,['comparison','comparison_type','method','platform','arm','context','proxy','cell_line','metric'],'PAIRED')
# Coverage/status counts include excluded and semantically unavailable rows.
counts=Counter(tuple(r[k] for k in ['task','platform','arm','context','proxy','cell_line','metric','status']) for r in ct)
write('CTSE_STATUS_COVERAGE.tsv',[dict(zip(['task','platform','arm','context','proxy','cell_line','metric','status','library_rows'],(*k,v))) for k,v in sorted(counts.items())])
dr=read('DRIFT_COMPOSITION_LONG.tsv');pr=read('REACTOME_COMPOSITION_LONG.tsv');ps=[]
dg=['task','method','platform','arm','context','cell_line','composition']
for key,rows in sorted(group(pr,dg).items()):
 vals=[abs(num(r['P'])) for r in rows if valid(r['P'])];ps.append(dict(zip(dg,key),metric='DIRECTIONAL_THRESHOLD_EXCEEDANCE',value=mean(v>.05 for v in vals),pathways_planned=len(rows),pathways_available=len(vals),median_abs_P=quantile(vals,.5),interpretation='DIAGNOSTIC_NOT_FPR_FDR_OR_SIGNIFICANCE'))
write('PROJECTION_COMPOSITION_SUMMARY.tsv',ps)
ds=[]
for rows,field in [(dr,'value'),(ps,'value')]:
 for key,rr in sorted(group(rows,['task','method','platform','arm','context','cell_line','metric']).items()):
  vv=[r[field] for r in rr];ds.append(dict(zip(['task','method','platform','arm','context','cell_line','metric'],key),planned=5,available=sum(valid(v) for v in vv),median=quantile(vv,.5),Q1=quantile(vv,.25),Q3=quantile(vv,.75),IQR=quantile(vv,.75)-quantile(vv,.25)))
write('DRIFT_AND_PROJECTION_SUMMARY.tsv',ds)
# Same support for relative A and P; no native cross-method subtraction.
dp=[]
for rr,field,metric,extra in [(dr,'value','DRIFT_AMPLITUDE_A',[]),(pr,'P','SIGNED_PROJECTION_P',['pathway'])]:
 keyfields=['platform','context','cell_line','composition']+extra
 bi={m:{tuple(r[k] for k in keyfields):r for r in rr if r['method']==m and not r['context'].endswith('_NATIVE')} for m in basenames}
 oi={tuple(r[k] for k in ['method']+keyfields):r for r in rr if r['arm']=='ORACLE_FRACTION_SENSITIVITY' and '_PAIR_' in r['context']}
 for a in rr:
  if a['method'] in basenames or a['context'].endswith('_NATIVE'):continue
  pairs=[(a['method']+'_MINUS_'+m,bi[m][tuple(a[k] for k in keyfields)]) for m in basenames]
  if a['arm']=='PRACTICAL_NONORACLE' and '_PAIR_' in a['context']:pairs.append((a['method']+'_PRACTICAL_MINUS_ORACLE',oi[tuple(a[k] for k in ['method']+keyfields)]))
  for comp,b in pairs:
   assert a['gene_count']==b['gene_count'];ok=valid(a[field]) and valid(b[field]);dp.append(dict(comparison=comp,method=a['method'],platform=a['platform'],arm=a['arm'],context=a['context'],cell_line=a['cell_line'],composition=a['composition'],pathway=a.get('pathway','NOT_APPLICABLE'),metric=metric,left=num(a[field]),right=num(b[field]),value=num(a[field])-num(b[field]) if ok else math.nan,status='AVAILABLE' if ok else 'PAIRED_ENDPOINT_UNAVAILABLE'))
write('DRIFT_PROJECTION_PAIRED_LONG.tsv',dp)
