import csv,pathlib,sys,os
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]/'repro'))
from archive import write_hashes,publish,sha256
R=pathlib.Path(sys.argv[1]);base=R/'audits'/'artificial_archive';base.mkdir()
src=base/'staging';src.mkdir();(src/'STATUS.tsv').write_text('task_id\tstatus\nartificial\tARTIFICIAL_ONLY\n');(src/'data.txt').write_text('completely artificial fixture; no CellBench/simulation values\n')
write_hashes(src);expected=sha256(src/'data.txt');dest=base/'published';publish(src,dest,{'STATUS.tsv','data.txt','hashes.sha256'})
checks=[('atomic_publish',dest.is_dir() and not src.exists()),('readonly_modes',all(not p.stat().st_mode & 0o222 for p in [dest,*dest.iterdir()])),('hash_preserved',sha256(dest/'data.txt')==expected)]
try:
 with open(dest/'data.txt','a') as f:f.write('SHOULD_FAIL')
 checks.append(('actual_write_blocked',False))
except PermissionError:checks.append(('actual_write_blocked',True))
duplicate=base/'duplicate';duplicate.mkdir();(duplicate/'STATUS.tsv').write_text('task_id\tstatus\nartificial\tARTIFICIAL_ONLY\n');write_hashes(duplicate)
try:publish(duplicate,dest,{'STATUS.tsv','hashes.sha256'});checks.append(('duplicate_blocked',False))
except FileExistsError:checks.append(('duplicate_blocked',True))
checks.append(('existing_result_not_overwritten',sha256(dest/'data.txt')==expected))
with open(R/'audits/ARCHIVE_MOCK.tsv','x',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['check','status']);w.writerows((k,'PASS' if v else 'FAIL') for k,v in checks)
assert all(v for _,v in checks)
print('ARCHIVE_MOCK_PASS',len(checks))
