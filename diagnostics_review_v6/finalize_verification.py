"""Verify report regeneration/publication protection after both replay checks."""
import sys,hashlib,json,shutil,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from run_benchmark import code_hashes


def main():
    source=ROOT/'results_v6';out=ROOT/'results_v6_validation'
    report=json.loads((source/'benchmark_verification_summary.json').read_text())
    assert report['status']=='completed' and report['source_sha256']==code_hashes()
    for code in ['benchmark_a','benchmark_b']:
        result=json.loads((out/(code+'_saved_verification.json')).read_text())
        assert result['status']=='saved_arrays_metrics_and_mc_precision_reproduced'
        assert result['source_sha256']==code_hashes()
    notes=source/'DIAGNOSTIC_EDITORIAL_NOTES.md'
    if not notes.exists():notes.write_text('Galat Monte Carlo menyatakan ketelitian estimasi dari campuran yang dipasang. '
        'Galat tersebut berbeda dari spread respons model, noise pengukuran, dan kalibrasi posterior.\n')
    notes_before=notes.read_bytes()
    command=[sys.executable,str(ROOT/'generate_pdf_manuscript.py'),'--results-dir',str(source),'--diagnostic-report']
    subprocess.run(command,check=True,cwd=ROOT)
    first=(source/'DIAGNOSTIC_VERIFICATION_REPORT.md').read_bytes()
    subprocess.run(command,check=True,cwd=ROOT)
    second=(source/'DIAGNOSTIC_VERIFICATION_REPORT.md').read_bytes()
    assert first==second and notes.read_bytes()==notes_before
    text=second.decode()
    assert notes_before.decode().strip() in text and 'Model tengah prior:' in text
    for code in ['benchmark_a','benchmark_b']:
        for warning in report[code]['baseline_comparison']['warnings']:assert warning in text
        assert str(report[code]['predictive_metrics']['misfit_mc_se']) in text
    notice=ROOT/'NASKAH_LENGKAP_PENELITIAN_PHD_MFJ_SMVI.md';notice_before=notice.read_bytes()
    blocked=subprocess.run([sys.executable,str(ROOT/'generate_pdf_manuscript.py'),'--results-dir',str(source),'--no-pdf'],capture_output=True,text=True,cwd=ROOT)
    assert blocked.returncode!=0 and 'Research manuscript blocked:' in blocked.stderr
    assert notice.read_bytes()==notice_before
    assert not (ROOT/'NASKAH_LENGKAP_PENELITIAN_PHD_MFJ_SMVI.pdf').exists()
    for name,sha in report['source_sha256'].items():
        assert hashlib.sha256((source/'source_snapshot'/name).read_bytes()).hexdigest()==sha
    record={'status':'passed','repeated_report_identical':True,'structured_prior_warnings_preserved':True,
        'separate_editorial_notes_preserved':True,'mc_precision_present':True,'research_publication_blocked':True,
        'original_withdrawal_notice_preserved':True,'source_snapshot_verified':True,
        'source_sha256':code_hashes(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (out/'report_and_publication_verification.json').write_text(json.dumps(record,indent=2))
    shutil.copy2('/tmp/mfj-v6-run.log',source/'execution.log')
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(source.iterdir()) if p.is_file()}
    (out/'verified_artifact_sha256.json').write_text(json.dumps(hashes,indent=2))
    print('Repeated report, critical warnings, notes, MC precision and publication protection: passed')

if __name__=='__main__':main()
