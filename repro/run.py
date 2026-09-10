"""Rebuild frozen evidence in a clean environment. Never fits a scientific model."""
import argparse,concurrent.futures,csv,json,os,subprocess,sys
from pathlib import Path
from verify import sha

COMPONENTS=['simulation','reference','cellbench','diagnostics']

def verify_inputs(data, component):
    rows=json.loads((data/'input-manifest.json').read_text())
    prefixes=[('simulation' if component=='reference' else component)+'/', 'expected/'+component+'/']
    selected=[r for r in rows if any(r['file'].startswith(p) for p in prefixes)
              or r['file']=='expected/expected-results.json']
    assert selected,'EMPTY_INPUT_INVENTORY'
    for r in selected:
        p=data/r['file']
        assert p.is_file() and not p.is_symlink() and sha(p)==r['sha256'],'INPUT_HASH_MISMATCH:'+r['file']
    return selected

def one(component, data, output):
    out=output/component
    assert not out.exists(),'OUTPUT_EXISTS_NO_OVERWRITE'
    before=verify_inputs(data,component)
    source=data/('simulation' if component=='reference' else component)
    subprocess.run(['Rscript','--vanilla',f'study/rebuild_{component}.R',str(source),str(out)],check=True)
    subprocess.run([sys.executable,f'study/summarize_{component}.py',str(out)],check=True)
    subprocess.run([sys.executable,'repro/verify.py',str(data),str(out),component],check=True)
    for r in before:assert sha(data/r['file'])==r['sha256'],'INPUT_CHANGED_DURING_REBUILD'
    return component

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--component',choices=['all']+COMPONENTS,default='all')
    p.add_argument('--jobs',type=int,default=1,choices=range(1,5))
    a=p.parse_args();data=a.data.resolve();output=a.output.resolve()
    assert (data/'input-manifest.json').is_file(),'FROZEN_DATA_MANIFEST_REQUIRED'
    assert not output.exists(),'OUTPUT_EXISTS_NO_OVERWRITE'
    output.mkdir(parents=True)
    env={k:os.environ.get(k,'') for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','PYTHONHASHSEED','TZ','LC_ALL']}
    (output/'RUNTIME.json').write_text(json.dumps(dict(python=sys.version,environment=env,model_execution='DISABLED',input_manifest_sha256=sha(data/'input-manifest.json')),indent=2)+'\n')
    components=COMPONENTS if a.component=='all' else [a.component]
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:
        done=list(pool.map(lambda c:one(c,data,output),components))
    (output/'REBUILD_COMPLETE.json').write_text(json.dumps(dict(status='PASS',components=done,new_model_fits=0),indent=2)+'\n')

if __name__=='__main__':main()
