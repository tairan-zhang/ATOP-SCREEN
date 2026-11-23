"""Exercise the packaged analysis worker with synthetic count data and its actual native tools."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

root = Path(__file__).resolve().parents[1]
engine = root / "packages/analysis-engine"
worker = engine / "inst/worker/analysis_worker.R"
if sys.platform == "darwin":
    runtime = root / "packaging/runtime/macos-arm64"
    command = [str(root / "scripts/packaged-rscript-macos-arm64.sh"), str(worker)]
    environment = dict(os.environ, ATOP_RUNTIME_ROOT=str(runtime),
                       ATOP_CPP_ENGINE_LIBRARY=str(runtime / "analysis/lib/atop_permutation_engine.so"))
else:
    command = ["Rscript", "--vanilla", str(worker)]
    environment = dict(os.environ)
environment["ATOP_ANALYSIS_ROOT"] = str(engine)

with tempfile.TemporaryDirectory() as folder:
    work = Path(folder)
    counts = work / "input matrix.tsv"
    with counts.open("w", newline="") as stream:
        table = csv.writer(stream, delimiter="\t")
        table.writerow(["guide","gene","T0","E1","E2","E3"])
        for i in range(90):
            table.writerow([f"g{i}", f"G{i//3}", 30+(3*i)%57, 12+(17*i)%101,
                            9+(13*i)%87, 10+(23*i)%123])
    cases = [
        ("paired", ["E1","E2"], ["T0","E3"], "cpp", False),
        ("auto", ["E1","E2","E3"], ["T0"], "cpp", True),
        ("group_means", ["E1"], ["T0","E3"], "cpp", False),
        ("group_means", ["E1","E2"], ["T0"], "r_parallel", False),
    ]
    requests = []
    for mode,endpoints,baselines,backend,mageck in cases:
        for seed in (17,17,85):
            idx = len(requests)
            params = dict(gRNA_col="guide",gene_col="gene",comparison_mode=mode,
                          condition1_replicate_cols=endpoints,condition2_replicate_cols=baselines,
                          N_perm=51,permutation_seed=seed,user_engine_choice=backend,run_mageck=mageck)
            requests.append(dict(id=str(idx),command="process_data",path=str(counts),parameters=params,
                                 artifactDir=str(work / str(idx))))
    proc = subprocess.run(command,cwd=engine,env=environment,text=True,capture_output=True,
                          input="".join(json.dumps(request)+"\n" for request in requests),check=True,timeout=600)
    replies = [row for line in proc.stdout.splitlines() if line.startswith("{") and "id" in (row:=json.loads(line))]
    assert len(replies)==len(requests),proc.stderr
    for idx,(request,reply) in enumerate(zip(requests,replies)):
        assert reply["id"]==str(idx) and reply["ok"],reply
        table=reply["result"]["tables"]
        settings=json.loads((work/str(idx)/"analysis-settings.json").read_text())
        assert settings["permutation_seed"]==request["parameters"]["permutation_seed"]
        backend=request["parameters"]["user_engine_choice"]
        assert settings["permutation_engine_used"]=={"cpp":"C++17","r_parallel":"R parallel"}[backend]
        if request["parameters"]["run_mageck"]:
            assert settings["comparison_mode"]=="shared_baseline" and len(table["mageckGenes"])==30
            assert len(table["mageckGuides"])==90 and (work/str(idx)/"mageck.gene_summary.txt").is_file()
        if idx%3:
            first=replies[idx-idx%3]["result"]["tables"]["genes"]
            assert (table["genes"]==first)==(idx%3==1)
print("Paired, shared baseline, group means, both seeded permutation engines and MAGeCK RRA passed.")
