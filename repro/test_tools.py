"""Artificial negative tests for the release gate and safe archive extraction."""
import csv,hashlib,io,json,sys,tarfile,tempfile
from pathlib import Path
from verify import compare,sha
from download_data import extract_archive
from privacy_guard import inspect

def main():
    tests=[]
    with tempfile.TemporaryDirectory() as temporary:
        root=Path(temporary);expected=root/'expected.tsv';actual=root/'actual.tsv'
        header=['task','status','count','value'];rows=[['a','OK','2','1.0'],['b','FAILED','0','NA']]
        def write(path,data):
            with path.open('w',newline='') as f:
                w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(header);w.writerows(data)
        write(expected,rows)
        record=dict(file='artificial.tsv',sha256=sha(expected),rows=2,columns=header,continuous_columns=['value'])
        def check(name,observed,want):
            write(actual,observed);status=compare(expected,actual,record)['status']
            assert status.startswith(want),(name,status,want);tests.append(name)
        check('exact',rows,'BYTE_IDENTICAL')
        check('row_order',rows[::-1],'EXACT_VALUES')
        check('declared_roundoff',[['a','OK','2','1.0000000000001'],rows[1]],'NUMERIC_MATCH')
        check('changed_result_blocks',[['a','OK','2','1.01'],rows[1]],'NUMERIC_VALUE_MISMATCH')
        check('status_change_blocks',[rows[0],['b','OK','0','NA']],'KEY_COUNT_STATUS')
        check('count_change_blocks',[['a','OK','3','1.0'],rows[1]],'KEY_COUNT_STATUS')
        check('failed_zero_imputation_blocks',[rows[0],['b','FAILED','0','0']],'MISSINGNESS_MASK')
        check('row_deletion_blocks',[rows[0]],'SCHEMA_OR_ROW_COUNT')
        check('duplicate_key_blocks',[rows[0],rows[0]],'DUPLICATE_SCIENTIFIC_KEYS')
        assert compare(expected,root/'missing.tsv',record)['status']=='MISSING';tests.append('missing_table_blocks')
        for kind,name in [('path','../escape'),('symlink','safe'),('unregistered','extra')]:
            archive=root/(kind+'.tar.gz');target=root/(kind+'-out');target.mkdir()
            with tarfile.open(archive,'w:gz') as f:
                info=tarfile.TarInfo(name)
                if kind=='symlink':info.type=tarfile.SYMTYPE;info.linkname='outside'
                else:info.size=1
                f.addfile(info,io.BytesIO(b'x') if info.isfile() else None)
            try:extract_archive(archive,target,{'safe'})
            except AssertionError:tests.append(kind+'_archive_blocked')
            else:raise AssertionError('UNSAFE_ARCHIVE_ACCEPTED')
        private=root/'privacy';private.mkdir()
        # Assemble artificial sensitive tokens; never include real credentials.
        (private/'example.txt').write_text('/'+'home'+'/artificial/private\n')
        assert inspect(private)['status']=='FAIL';tests.append('private_path_detection')
        (private/'example.txt').write_text('C'+':/'+'Users'+'/artificial/private\n')
        assert inspect(private)['status']=='FAIL';tests.append('windows_private_path_detection')
        (private/'example.txt').write_text('"$PWD:/work:ro"\n')
        assert inspect(private)['status']=='PASS';tests.append('public_container_mount_not_personal_drive')
    print(json.dumps(dict(status='PASS',artificial_only=True,checks=tests,model_fits=0),indent=2))

if __name__=='__main__':main()
