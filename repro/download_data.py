"""Download/hash-check versioned scientific assets and extract without overwrite."""
import argparse,hashlib,json,os,tarfile,urllib.request
from pathlib import Path,PurePosixPath

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def extract_archive(path,root,allowed):
    seen=[]
    with tarfile.open(path,'r|gz') as archive:
        for member in archive:
            name=PurePosixPath(member.name)
            assert member.isfile() and not name.is_absolute() and '..' not in name.parts,'UNSAFE_ARCHIVE_MEMBER'
            assert member.name in allowed and member.name not in seen,'UNREGISTERED_OR_DUPLICATE_MEMBER'
            seen.append(member.name);target=root.joinpath(*name.parts)
            target.parent.mkdir(parents=True,exist_ok=True)
            assert not target.exists() and not target.is_symlink(),'NO_OVERWRITE'
            with archive.extractfile(member) as inp,target.open('xb') as out:
                for b in iter(lambda:inp.read(1048576),b''):out.write(b)
    assert set(seen)==set(allowed),'ARCHIVE_INVENTORY_MISMATCH'

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--destination',type=Path,required=True)
    p.add_argument('--asset-directory',type=Path,help='Optional existing verified archive directory; disables downloads')
    p.add_argument('--cache',type=Path,default=Path('.asset-cache'))
    p.add_argument('--manifest',type=Path,default=Path(__file__).with_name('data-assets.json'))
    a=p.parse_args();spec=json.loads(a.manifest.read_text());dest=a.destination.resolve()
    assert not dest.exists() and not dest.is_symlink(),'DESTINATION_EXISTS_NO_OVERWRITE'
    stage=dest.with_name(dest.name+'.building')
    assert not stage.exists(),'STAGING_EXISTS_REQUIRES_INSPECTION'
    cache=a.asset_directory if a.asset_directory else a.cache
    if not a.asset_directory:cache.mkdir(parents=True,exist_ok=True)
    archives=[]
    for row in spec['assets']:
        path=cache/row['name'];assert path.parent==cache,'UNSAFE_ASSET_NAME'
        if not path.exists():
            assert not a.asset_directory,'OFFLINE_ASSET_MISSING'
            temporary=path.with_name(path.name+'.part');assert not temporary.exists(),'PARTIAL_DOWNLOAD_REQUIRES_INSPECTION'
            url='https://github.com/'+spec['repository']+'/releases/download/'+spec['tag']+'/'+row['name']
            with urllib.request.urlopen(url,timeout=180) as r,temporary.open('xb') as f:
                for b in iter(lambda:r.read(1048576),b''):f.write(b)
            assert sha(temporary)==row['sha256'],'DOWNLOAD_HASH_MISMATCH'
            temporary.rename(path)
        assert not path.is_symlink() and path.stat().st_size==row['bytes'] and sha(path)==row['sha256'],'ASSET_HASH_MISMATCH'
        archives.append((path,row));print('ASSET_VERIFIED',row['name'],flush=True)
    stage.mkdir(parents=True)
    for path,row in archives:extract_archive(path,stage,row['files'])
    manifest=stage/'input-manifest.json'
    assert sha(manifest)==spec['data_manifest_sha256'],'DATA_MANIFEST_HASH_MISMATCH'
    records=json.loads(manifest.read_text())
    files={x.relative_to(stage).as_posix() for x in stage.rglob('*') if x.is_file()}
    assert files=={r['file'] for r in records}|{'input-manifest.json'},'UNCLOSED_DATA_INVENTORY'
    for row in records:assert sha(stage/row['file'])==row['sha256'],'EXTRACTED_INPUT_HASH_MISMATCH'
    # Kernel-enforced no-replace publication; do not race an existing directory.
    from archive import rename_noreplace
    rename_noreplace(stage,dest)
    for path in dest.rglob('*'):path.chmod(0o555 if path.is_dir() else 0o444)
    dest.chmod(0o555)
    print('FROZEN_DATA_READY',len(files),'files')

if __name__=='__main__':main()
