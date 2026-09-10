"""Exact frozen bounded diagnostic summaries; no sampler execution."""
import collections,csv,pathlib,sys,math
A=pathlib.Path(sys.argv[1]);out=A/'SUMMARIES';out.mkdir(exist_ok=False)
def read(p):
 with open(p) as f:yield from csv.DictReader(f,delimiter='\t')
def write(p,rows,fields=None):
 rows=list(rows)
 with open(p,'x',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def num(v):
 try:return float(v)
 except (ValueError,TypeError):return float('nan')
def f(v):return format(v,'.17g') if math.isfinite(v) else 'NA'
inventory=list(read(A/'CHAIN_EXECUTION_INVENTORY.tsv'));assert len(inventory)==1920
chain_counts=collections.Counter(x['status'] for x in inventory)
fits=sum(int(x['fit_count']) for x in inventory if x['fit_count'] not in ['NA',''])
grouped={}
for x in read(A/'ALL_PARAMETER_DIAGNOSTICS.tsv'):
 key=tuple(x[k] for k in ['unit_id','stage','schedule','quantity'])
 g=grouped.setdefault(key,dict(n=0,passing=0,genes=set(),statuses=collections.Counter(),rhats=[],bulk=[],tail=[]))
 g['n']+=1;g['passing']+=x['status']=='PASS';g['genes'].add(x['gene_id']);g['statuses'][x['status']]+=1
 for dest,key2 in [('rhats','rhat'),('bulk','ess_bulk'),('tail','ess_tail')]:
  v=num(x[key2])
  if math.isfinite(v):g[dest].append(v)
summary=[]
for key,g in sorted(grouped.items()):
 summary.append(dict(zip(['unit_id','stage','schedule','quantity'],key))|dict(
  gene_count=len(g['genes']),parameter_records=g['n'],passing_records=g['passing'],
  diagnostic_pass_fraction=f(g['passing']/g['n']),max_rhat=f(max(g['rhats'],default=float('nan'))),
  min_bulk_ess=f(min(g['bulk'],default=float('nan'))),min_tail_ess=f(min(g['tail'],default=float('nan'))),
  statuses=';'.join(k+':'+str(v) for k,v in sorted(g['statuses'].items())),
  interpretation='DESCRIPTIVE_PARAMETER_SCREEN_NOT_INDEPENDENT_REPLICATES'))
write(out/'NUMERICAL_DIAGNOSTIC_COVERAGE.tsv',summary)
paired={}
for x in read(A/'POSTERIOR_SUMMARIES_LONG.tsv'):
 key=tuple(x[k] for k in ['unit_id','stage','gene_id','quantity','parameter'])
 d=paired.setdefault(key,{});assert x['schedule'] not in d;d[x['schedule']]=x
comparison=[]
for key,d in sorted(paired.items()):
 a,b=d.get('DEFAULT',{}),d.get('LONG',{})
 z=dict(zip(['unit_id','stage','gene_id','quantity','parameter'],key))
 for tag,x in [('default',a),('long',b)]:
  for name in ['posterior_mean','posterior_sd','mcse_mean','status']:z[tag+'_'+name]=x.get(name,'MISSING')
 z['long_minus_default_mean']=f(num(b.get('posterior_mean'))-num(a.get('posterior_mean')))
 z['interpretation']='FIXED_POSTERIOR_NUMERICAL_SENSITIVITY_NOT_PERFORMANCE'
 comparison.append(z)
write(out/'POSTERIOR_LENGTH_COMPARISON.tsv',comparison)

print('DIAGNOSTIC_SUMMARY_REBUILT_NO_SAMPLING')
