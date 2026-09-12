"""Fail-closed checks against the predeclared frozen scientific table inventory."""
import argparse,csv,hashlib,json,math
from pathlib import Path
import numpy as np
import pandas as pd

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def clean_text(frame):
    return frame.replace({'TRUE':'True','FALSE':'False','true':'True','false':'False'})

def compare(expected, observed, record):
    result=dict(file=record['file'],expected_sha256=record['sha256'],observed_sha256='',
                status='MISSING',rows_expected=record['rows'],rows_observed=0,
                differing_numeric_cells=0,max_absolute_difference=0.0)
    assert sha(expected)==record['sha256'],'EXPECTED_ASSET_HASH_DRIFT'
    if not observed.is_file():return result
    result['observed_sha256']=sha(observed)
    if result['observed_sha256']==record['sha256']:
        result.update(status='BYTE_IDENTICAL',rows_observed=record['rows']);return result
    e=clean_text(pd.read_csv(expected,sep='\t',dtype=str,keep_default_na=False))
    a=clean_text(pd.read_csv(observed,sep='\t',dtype=str,keep_default_na=False))
    result['rows_observed']=len(a)
    if set(a.columns)!=set(e.columns) or len(a)!=len(e):
        result['status']='SCHEMA_OR_ROW_COUNT_MISMATCH';return result
    nums=record['continuous_columns'];keys=[c for c in e.columns if c not in nums]
    # Every non-continuous field, including counts, statuses and support IDs,
    # participates in reconciliation and must match exactly.
    if e.duplicated(keys).any() or a.duplicated(keys).any():
        result['status']='DUPLICATE_SCIENTIFIC_KEYS';return result
    e=e.sort_values(keys,kind='stable').reset_index(drop=True)
    a=a.sort_values(keys,kind='stable').reset_index(drop=True)
    if not e[keys].equals(a[keys]):
        result['status']='KEY_COUNT_STATUS_OR_ELIGIBILITY_MISMATCH';return result
    exact=True
    for c in nums:
        missing_e=e[c].isin(['NA','NaN','nan','']);missing_a=a[c].isin(['NA','NaN','nan',''])
        if not missing_e.equals(missing_a):
            result['status']='MISSINGNESS_MASK_MISMATCH:'+c;return result
        ev=pd.to_numeric(e.loc[~missing_e,c],errors='raise').to_numpy(dtype=float)
        av=pd.to_numeric(a.loc[~missing_a,c],errors='raise').to_numpy(dtype=float)
        if not (np.isfinite(ev).all() and np.isfinite(av).all()):
            result['status']='NONFINITE_UNDECLARED_VALUE:'+c;return result
        diff=np.abs(av-ev);n=int(np.count_nonzero(diff));exact &= n==0
        result['differing_numeric_cells']+=n
        result['max_absolute_difference']=max(result['max_absolute_difference'],float(diff.max(initial=0)))
        if not np.all(diff<=1e-12+1e-12*np.abs(ev)):
            result['status']='NUMERIC_VALUE_MISMATCH:'+c;return result
    result['status']='EXACT_VALUES_ROW_ORDER_OR_SERIALIZATION_ONLY' if exact else 'NUMERIC_MATCH_WITH_DECLARED_ROUNDOFF'
    return result

def main(data,output,component):
    records=json.loads((data/'expected/expected-results.json').read_text())
    selected=[r for r in records if r['component']==component];assert selected
    results=[compare(data/'expected'/component/r['file'],output/r['file'],r) for r in selected]
    passed={'BYTE_IDENTICAL','EXACT_VALUES_ROW_ORDER_OR_SERIALIZATION_ONLY','NUMERIC_MATCH_WITH_DECLARED_ROUNDOFF'}
    status='PASS' if all(r['status'] in passed for r in results) else 'FAIL'
    with (output/'REPRODUCTION_AUDIT.tsv').open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(results)
    report=dict(component=component,status=status,tables_expected=len(selected),tables_checked=len(results),
                new_model_fits=0,policy='repro/COMPARISON_POLICY.md',results=results)
    (output/'REPRODUCTION_STATUS.json').write_text(json.dumps(report,indent=2)+'\n')
    print(component,status,len(results),'tables',flush=True)
    if status!='PASS':
        for r in results:
            if r['status'] not in passed:print(r['file'],r['status'],flush=True)
        raise SystemExit(1)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('data',type=Path);p.add_argument('output',type=Path);p.add_argument('component')
    a=p.parse_args();main(a.data,a.output,a.component)
