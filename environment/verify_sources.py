import hashlib,json
from pathlib import Path
root=Path('/tmp')
rows=json.loads((root/'r-source-lock.json').read_text())
assert {p.name for p in (root/'r-sources').iterdir()}=={Path(r['file']).name for r in rows}
for r in rows:
    p=root/r['file']
    assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],'R_SOURCE_HASH_DRIFT:'+r['file']
print('All pinned R source hashes matched.')
