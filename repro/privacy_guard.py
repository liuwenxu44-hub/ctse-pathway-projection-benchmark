"""Reject private paths, credentials and operational files in a public tree."""
import argparse,hashlib,json,re
from pathlib import Path

PATTERNS={
 'personal_or_server_path':re.compile(r'/(?:home|Users|mnt)/|/data/[A-Za-z][A-Za-z0-9_-]*/|(?<![A-Za-z0-9_])[A-Z]:[\\/](?:Users|work|proj)',re.I),
 'private_IP':re.compile(r'\b(?:10\.[0-9]{1,3}|192\.168|172\.(?:1[6-9]|2[0-9]|3[01]))\.[0-9]{1,3}\.[0-9]{1,3}\b'),
 'github_credential':re.compile(r'\b(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,})\b'),
 'private_key':re.compile(r'-----BEGIN [A-Z ]*PRIVATE KEY-----'),
 'embedded_auth_url':re.compile(r'https?://[^/\s]+:[^/\s]+@'),
}
FORBIDDEN={'.env','.Rhistory','.RData','id_rsa','id_ed25519','config.toml'}
def inspect(root):
    errors=[];files=0
    for p in sorted(root.rglob('*')):
        rel=p.relative_to(root)
        if '.git' in rel.parts or '__pycache__' in rel.parts:continue
        if p.is_symlink():errors.append(dict(file=str(rel),rule='SYMLINK_NOT_ALLOWED'));continue
        if not p.is_file():continue
        files+=1
        if p.name in FORBIDDEN or p.suffix.lower() in {'.log','.docx','.pdf','.tif','.tiff','.ai','.psd'}:
            errors.append(dict(file=str(rel),rule='OPERATIONAL_OR_MANUSCRIPT_ARTIFACT_NOT_IN_RELEASE'))
        if p.suffix.lower() not in {'.r','.py','.md','.tsv','.csv','.json','.yml','.yaml','.txt'} and p.name not in {'Dockerfile','LICENSE','.gitignore'}:continue
        try:content=p.read_text(encoding='utf-8')
        except UnicodeDecodeError:errors.append(dict(file=str(rel),rule='TEXT_NOT_UTF8'));continue
        for name,pattern in PATTERNS.items():
            if pattern.search(content):errors.append(dict(file=str(rel),rule=name))
    return dict(status='PASS' if not errors else 'FAIL',files=files,findings=errors)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('root',type=Path);a.add_argument('--report',type=Path);args=a.parse_args()
    r=inspect(args.root)
    if args.report:
        with args.report.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(r,indent=2));raise SystemExit(r['status']!='PASS')
