#!/usr/bin/env bash
# Fifth fit (setup in fit_config.json): the stages one after another; a failed stage stops the chain.
# Skipped on purpose, as in the third and fourth fits: 'refit' (repaired a solver setting of the first run) and 'throughput'.
# Writes provenance.txt and code_vs_commit.patch when it starts. Run it inside tmux: tmux new -s fit5 bash <this file>
cd /home/yeray/projects/atom-to-reactor || exit 1
PY=/home/yeray/miniforge3/envs/atom-to-reactor/bin/python
OUT=notebooks/results/05_fit5
DATA=$(grep '^ATOM_DATA_DIR=' .env | cut -d= -f2-)
if [ ! -f "$OUT/provenance.txt" ]; then
    git diff HEAD -- kinetics scripts tests data > "$OUT/code_vs_commit.patch"
    {
        echo "Fifth fit (setup: $OUT/fit_config.json), launched $(date '+%Y-%m-%d %H:%M:%S %Z') on $(hostname)"
        echo "code: branch $(git rev-parse --abbrev-ref HEAD), commit $(git rev-parse HEAD) plus the uncommitted changes in code_vs_commit.patch (kinetics/, scripts/, tests/, data/)"
        echo "untracked code files: $(git ls-files --others --exclude-standard -- kinetics scripts tests | tr '\n' ' ')"
        echo "observables: $(md5sum "$DATA/data/lab_observables.json")"
        echo "snapshot (solvation): $(md5sum "$DATA/data/tank_api_snapshot_2026-10-08.json")"
        echo "snapshot (gas phase, geometries, shieldings; identical in the two files): $(md5sum "$DATA/data/tank_api_snapshot_2026-09-24.json")"
        echo "stages: registered recovery baselines extended night traces predictions (14 workers); refit and throughput skipped"
    } > "$OUT/provenance.txt"
    # the untracked code files are copied so the run can be rebuilt
    mkdir -p "$OUT/code_untracked"
    git ls-files --others --exclude-standard -- kinetics scripts tests | while read -r f; do
        mkdir -p "$OUT/code_untracked/$(dirname "$f")"; cp "$f" "$OUT/code_untracked/$f"; done
else
    stamp=$(date +%y%m%d_%H%M)
    git diff HEAD -- kinetics scripts tests data > "$OUT/code_vs_commit_restart_$stamp.patch"
    cp "$OUT/fit_config.json" "$OUT/fit_config_restart_$stamp.json"
    echo "[$(date +%T)] chain restarted at commit $(git rev-parse HEAD); finished jobs are read from their files; code and setup at the restart: code_vs_commit_restart_$stamp.patch, fit_config_restart_$stamp.json" >> "$OUT/provenance.txt"
fi
for stage in registered recovery baselines extended night traces predictions; do
    "$PY" scripts/run_fit.py "$stage" --workers 14 --results "$OUT" >> "$OUT/fit5.out" 2>&1
    code=$?
    if [ $code -ne 0 ]; then
        echo "[$(date +%T)] CHAIN STOPPED: stage $stage failed (exit $code), see $OUT/fit5.out" | tee -a "$OUT/fit5.out" "$OUT/run.log"
        exit 1
    fi
done
echo "[$(date +%T)] CHAIN FINISHED: all stages done" | tee -a "$OUT/fit5.out" "$OUT/run.log"
