#!/usr/bin/env bash
# Eighth fit (setup in fit_config.json, 2026-10-23): M1, M5 and M5-first, each with Marcus and with P. Broqvist's
# rule ('+BEP'), at 10 min and with free ages between 10 min and 7 days: 12 fits, their recovery and profiles.
# Written to fill every core of the cloud VM (cloud.sh): every independent job is its own process.
#
#   1 registered                 the registered models and the computed energies (seconds)
#   2 extended                   the two readings side by side, one process each (inside a reading the structures
#                                go one after the other: each starts from the smaller ones)
#   3 recovery and night         one process per structure and reading for each of the two, all side by side
#                                (24 processes): the recovery (the fit as truth, 2 seeds) and the profiles
#   4 extended, recovery, night  again without limits: every result is read back, the tables are written
#   5 traces, predictions        the trace rule and the checks on data not fitted
#
#   tmux new -s fit8 bash notebooks/results/05_fit8/run_chain.sh          the run
#   tmux new -s fit8 bash notebooks/results/05_fit8/run_chain.sh smoke    the same chain with --quick in 05_fit8_smoke
#
# Every job writes its own file and a job already on disk is read back, not redone: the chain can be stopped and
# started again at any point (only the jobs in progress are lost). A failed process stops the chain after its
# step. Writes provenance.txt and code_vs_commit.patch when it starts.
#
# CODE is the worktree of branch fit/ten-min-free (the seventh fit's code plus the free-age range and --readings).
CODE=/home/yeray/projects/atom-to-reactor-fit8
OUT=/home/yeray/projects/atom-to-reactor/notebooks/results/05_fit8
PY=/home/yeray/miniforge3/envs/atom-to-reactor/bin/python
CPUS=$(nproc)
STRUCTURES=("M1" "M1+BEP" "M5" "M5+BEP" "M5-first" "M5-first+BEP")
READINGS=("10 min" "free")
QUICK=()
if [ "$1" = "smoke" ]; then
    mkdir -p "${OUT}_smoke"
    cp "$OUT/fit_config.json" "${OUT}_smoke/"
    OUT=${OUT}_smoke
    QUICK=(--quick)
fi
cd "$CODE" || exit 1
DATA=$(grep '^ATOM_DATA_DIR=' .env | cut -d= -f2-)
if [ ! -f "$OUT/provenance.txt" ]; then
    git diff HEAD -- kinetics scripts tests data > "$OUT/code_vs_commit.patch"
    {
        echo "Eighth fit (setup: $OUT/fit_config.json), launched $(date '+%Y-%m-%d %H:%M:%S %Z') on $(hostname), $CPUS CPUs"
        echo "code: $CODE, branch $(git rev-parse --abbrev-ref HEAD), commit $(git rev-parse HEAD) plus the uncommitted changes in code_vs_commit.patch (kinetics/, scripts/, tests/, data/)"
        echo "untracked code files: $(git ls-files --others --exclude-standard -- kinetics scripts tests | tr '\n' ' ')"
        echo "observables: $(md5sum "$DATA/data/lab_observables.json")"
        echo "snapshot (solvation): $(md5sum "$DATA/data/tank_api_snapshot_2026-10-08.json")"
        echo "snapshot (gas phase, geometries, shieldings; identical in the two files): $(md5sum "$DATA/data/tank_api_snapshot_2026-09-24.json")"
        echo "starts from: the fits of 05_fit7 (fit_config.json start_from), as they were at the launch"
        echo "steps: registered; extended per reading; recovery and night per structure and reading; extended, recovery and night without limits; traces; predictions${QUICK:+ (smoke: --quick)}"
    } > "$OUT/provenance.txt"
    mkdir -p "$OUT/code_untracked"
    git ls-files --others --exclude-standard -- kinetics scripts tests | while read -r f; do
        mkdir -p "$OUT/code_untracked/$(dirname "$f")"; cp "$f" "$OUT/code_untracked/$f"; done
else
    stamp=$(date +%y%m%d_%H%M)
    git diff HEAD -- kinetics scripts tests data > "$OUT/code_vs_commit_restart_$stamp.patch"
    cp "$OUT/fit_config.json" "$OUT/fit_config_restart_$stamp.json"
    echo "[$(date +%T)] chain restarted on $(hostname) ($CPUS CPUs) at commit $(git rev-parse HEAD); finished jobs are read from their files; code and setup at the restart: code_vs_commit_restart_$stamp.patch, fit_config_restart_$stamp.json" >> "$OUT/provenance.txt"
fi

PIDS=()
start() {            # one run_fit.py call in the background; its output to its own file under logs/
    local name=$1; shift
    mkdir -p "$OUT/logs"
    "$PY" scripts/run_fit.py "$@" "${QUICK[@]}" --results "$OUT" > "$OUT/logs/$name.out" 2>&1 &
    PIDS+=($!)
}
finish() {           # waits for every call started since the last finish; stops the chain if one failed
    local failed=0
    for pid in "${PIDS[@]}"; do wait "$pid" || failed=1; done
    PIDS=()
    if [ $failed -ne 0 ]; then
        echo "[$(date +%T)] CHAIN STOPPED at $1: a process failed, see $OUT/logs/" | tee -a "$OUT/run.log"
        exit 1
    fi
    echo "[$(date +%T)] STEP DONE: $1" | tee -a "$OUT/run.log"
}
tag() { echo "$1" | tr ' +' '_p'; }

start registered registered --workers "$CPUS"
finish "registered"
for reading in "${READINGS[@]}"; do
    start "extended_$(tag "$reading")" extended --readings "$reading" --workers $((CPUS / ${#READINGS[@]}))
done
finish "extended, per reading"
# 24 processes: each can use as many workers as it has tasks, up to a third of the cores (profile chains, starts)
for reading in "${READINGS[@]}"; do
    for structure in "${STRUCTURES[@]}"; do
        start "recovery_$(tag "$structure")_$(tag "$reading")" recovery --structures "$structure" --readings "$reading" --workers $((CPUS / 3))
        start "night_$(tag "$structure")_$(tag "$reading")" night --structures "$structure" --readings "$reading" --workers $((CPUS / 3))
    done
done
finish "recovery and night, per structure and reading"
for stage in extended recovery night; do
    start "$stage" "$stage" --workers "$CPUS"
    finish "$stage, every structure and reading (tables)"
done
for stage in traces predictions; do
    start "$stage" "$stage" --workers "$CPUS"
    finish "$stage"
done
echo "[$(date +%T)] CHAIN FINISHED" | tee -a "$OUT/run.log"
