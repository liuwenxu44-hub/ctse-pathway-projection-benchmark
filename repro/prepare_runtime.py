"""Prepare the pinned Docker build context; no host installation or model fitting."""
import argparse,csv,json,subprocess,sys,urllib.request
from pathlib import Path
from fetch_environment import sha

def download(url,path,expected):
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        assert not path.is_symlink() and sha(path)==expected,'EXISTING_RUNTIME_FILE_HASH_MISMATCH'
        return
    temporary=path.with_suffix(path.suffix+'.part')
    assert not temporary.exists(),'PARTIAL_DOWNLOAD_REQUIRES_INSPECTION'
    with urllib.request.urlopen(url,timeout=120) as r,temporary.open('xb') as f:
        while block:=r.read(1048576):f.write(block)
    assert sha(temporary)==expected,'RUNTIME_DOWNLOAD_HASH_MISMATCH'
    temporary.rename(path)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('environment',type=Path)
    env=parser.parse_args().environment.resolve()
    subprocess.run([sys.executable,str(Path(__file__).with_name('fetch_environment.py')),str(env)],check=True)
    download('https://cloud.r-project.org/bin/linux/ubuntu/noble-cran40/r-base-core_4.5.3-1.2404.0_amd64.deb',
        env/'debs/r-base-core_4.5.3-1.2404.0_amd64.deb','4066b325cc0e3e22ad7c786576fb8bef555ea0c09bd58d7a3f283969d87b9aa7')
    for row in csv.DictReader((env/'python-wheels.tsv').open(),delimiter='\t'):
        path=env/'wheels'/row['filename']
        if path.exists():
            assert not path.is_symlink() and sha(path)==row['sha256'],'EXISTING_WHEEL_HASH_MISMATCH'
            continue
        with urllib.request.urlopen('https://pypi.org/pypi/'+row['package']+'/'+row['version']+'/json',timeout=60) as r:
            info=json.load(r)
        matches=[x for x in info['urls'] if x['filename']==row['filename'] and x['digests']['sha256']==row['sha256']]
        assert len(matches)==1,'PINNED_OFFICIAL_WHEEL_NOT_FOUND'
        download(matches[0]['url'],path,row['sha256'])
    print('RUNTIME_CONTEXT_HASH_VERIFIED_NO_HOST_INSTALL')

if __name__=='__main__':main()
