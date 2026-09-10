"""Download exact third-party runtime sources, without installing on the host."""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import urllib.error
import urllib.request
from pathlib import Path

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def download(row, dest):
    name=f"{row['package']}_{row['version']}.tar.gz"
    p=dest/name
    urls=[f"https://cloud.r-project.org/src/contrib/Archive/{row['package']}/{name}",
          f"https://cloud.r-project.org/src/contrib/{name}"]
    for url in urls:
        try:
            if not p.exists():
                with urllib.request.urlopen(url,timeout=90) as r, p.with_suffix('.part').open('xb') as f:
                    while b:=r.read(1024*1024):f.write(b)
                p.with_suffix('.part').rename(p)
            return dict(package=row['package'],version=row['version'],license=row['license'],
                        file='r-sources/'+name,url=url,sha256=sha(p))
        except urllib.error.HTTPError as e:
            if e.code!=404:raise
    raise RuntimeError('Exact CRAN package unavailable: '+name)

def main():
    a=argparse.ArgumentParser();a.add_argument('environment',type=Path);args=a.parse_args()
    env=args.environment;dest=env/'r-sources';dest.mkdir(exist_ok=True)
    rows=list(csv.DictReader((env/'r-packages.tsv').open(),delimiter='\t'))
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        records=list(pool.map(lambda r:download(r,dest),rows))
    target=env/'r-source-lock.json'
    if target.exists():
        old=json.loads(target.read_text())
        assert {x['file']:x['sha256'] for x in old}=={x['file']:x['sha256'] for x in records},'SOURCE_HASH_DRIFT'
    else:target.write_text(json.dumps(records,indent=2)+'\n')
    print('Exact R sources:',len(records))

if __name__=='__main__':main()
