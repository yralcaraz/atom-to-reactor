#!/usr/bin/env bash
# Seventh fit (setup in fit_config.json; MET 2026-10-22 "Simple M1 fit at near-zero ages"): M1, M5 and M5-first,
# each with Marcus and with P. Broqvist's barrier rule ('+BEP'), ages of 0, 5 and 10 min: 18 fits, each analysed.
#
# Phase 1, M1 and M1+BEP (to finish within 10 h): every stage limited to them with --structures
#   registered    the registered models and the computed energies (seconds)
#   extended      the 6 fits
#   night         recovery (the fit as truth, 2 seeds), profiles in every reading, leave-one-out, 8 sensitivities,
#                 solver tolerance at 0 min
#   traces        a trace of TMSOH at mixing: each fit unchanged, and fitted again where it fits
#   predictions   tube B, the late 13C of tube A, TMSPa alone from its first spectrum, paper observations,
#                 storage at 25 °C, bench protocol, trajectories
# Phase 2, M5, M5-first and their '+BEP' forms (no time limit): the same stages without --structures; the M1 results
# are read back from their files, and the summary tables are written again with every structure.
#
#   tmux new -s fit7 bash notebooks/results/05_fit7/run_chain.sh       both phases
#   tmux new -s fit7 bash notebooks/results/05_fit7/run_chain.sh M1    phase 1 only; run it again without M1 to go on
#
# Every job writes its own file and a job already on disk is read back, not redone: the chain can be stopped and
# started again at any point. A failed stage stops the chain. Writes provenance.txt and code_vs_commit.patch when
# it starts.
#
# CODE is the worktree of branch fit/near-zero-ages, where the code of this fit lives (near-zero readings, M5,
# '+BEP', --structures, the night generalised to several structures); the main repository keeps the code of the
# sixth fit. Once the branch is merged, set CODE to /home/yeray/projects/atom-to-reactor.
CODE=/home/yeray/projects/atom-to-reactor-fit7
OUT=/home/yeray/projects/atom-to-reactor/notebooks/results/05_fit7
PY=/home/yeray/miniforge3/envs/atom-to-reactor/bin/python
WORKERS=14
M1="M1,M1+BEP"
cd "$CODE" || exit 1
DATA=$(grep '^ATOM_DATA_DIR=' .env | cut -d= -f2-)
if [ ! -f "$OUT/provenance.txt" ]; then
    git diff HEAD -- kinetics scripts tests data > "$OUT/code_vs_commit.patch"
    {
        echo "Seventh fit (setup: $OUT/fit_config.json), launched $(date '+%Y-%m-%d %H:%M:%S %Z') on $(hostname)"
        echo "code: $CODE, branch $(git rev-parse --abbrev-ref HEAD), commit $(git rev-parse HEAD) plus the uncommitted changes in code_vs_commit.patch (kinetics/, scripts/, tests/, data/)"
        echo "untracked code files: $(git ls-files --others --exclude-standard -- kinetics scripts tests | tr '\n' ' ')"
        echo "observables: $(md5sum "$DATA/data/lab_observables.json")"
        echo "snapshot (solvation): $(md5sum "$DATA/data/tank_api_snapshot_2026-10-08.json")"
        echo "snapshot (gas phase, geometries, shieldings; identical in the two files): $(md5sum "$DATA/data/tank_api_snapshot_2026-09-24.json")"
        echo "stages: phase 1 ($M1) registered extended night traces predictions; phase 2 (every structure) extended night traces predictions ($WORKERS workers)"
    } > "$OUT/provenance.txt"
    # the untracked code files are copied so the run can be rebuilt
    mkdir -p "$OUT/code_untracked"
    git ls-files --others --exclude-standard -- kinetics scripts tests | while read -r f; do
        mkdir -p "$OUT/code_untracked/$(dirname "$f")"; cp "$f" "$OUT/code_untracked/$f"; done
else
    stamp=$(date +%y%m%d_%H%M)
    git diff HEAD -- kinetics scripts tests data > "$OUT/code_vs_commit_restart_$stamp.patch"
    cp "$OUT/fit_config.json" "$OUT/fit_config_restart_$stamp.json"
    echo "[$(date +%T)] chain restarted at commit $(git rev-parse HEAD) (argument: ${1:-none}); finished jobs are read from their files; code and setup at the restart: code_vs_commit_restart_$stamp.patch, fit_config_restart_$stamp.json" >> "$OUT/provenance.txt"
fi

run() {
    "$PY" scripts/run_fit.py "$@" --workers "$WORKERS" --results "$OUT" >> "$OUT/fit7.out" 2>&1
    code=$?
    if [ $code -ne 0 ]; then
        echo "[$(date +%T)] CHAIN STOPPED: run_fit.py $* failed (exit $code), see $OUT/fit7.out" | tee -a "$OUT/fit7.out" "$OUT/run.log"
        exit 1
    fi
}

# Phase 1: M1 and M1+BEP
run registered
for stage in extended night traces predictions; do
    run "$stage" --structures "$M1"
done
echo "[$(date +%T)] PHASE 1 DONE (M1 and M1+BEP): fits, profiles, leave-one-out, sensitivities, traces and predictions" | tee -a "$OUT/fit7.out" "$OUT/run.log"
if [ "$1" = "M1" ]; then
    exit 0
fi
# Phase 2: M5, M5-first and their '+BEP' forms, every stage; M1 is read back
for stage in extended night traces predictions; do
    run "$stage"
done
echo "[$(date +%T)] CHAIN FINISHED: both phases done" | tee -a "$OUT/fit7.out" "$OUT/run.log"
