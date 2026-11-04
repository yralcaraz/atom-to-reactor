#!/usr/bin/env bash
# Eighth fit, M1 tube by tube (setup in fit_config.json; P. Broqvist 2026-10-23: one experiment at a time): the steps
# S1 glovebox, S2 tube A, S3 TMSPa alone, S2b tube A again, S4 0.5 % H2O, S5 2 % H2O, at 10 min, 1 h and 8 h; the
# free barriers of every step profiled; after every step every tube scored, the 2 % tube also as P. Broqvist reads it.
#
#   tmux new -s fit8 bash notebooks/results/fit8_M1-steps/run_chain.sh          the run
#   tmux new -s fit8 bash notebooks/results/fit8_M1-steps/run_chain.sh smoke    the same with --quick in fit8_M1-steps_smoke
#
# registered writes the computed energies (computed_energies.csv); steps does the rest. Every job writes its own
# file and a job already on disk is read back, not redone: the chain can be stopped and started again. Writes
# provenance.txt and code_vs_commit.patch when it starts.
#
# CODE is the worktree of branch fit/m1-steps (the seventh fit's code plus the ages 1 h and 8 h and the steps stage);
# the 05_fit8 setup (free ages, branch fit/ten-min-free) is another run. Once the branch is merged, set CODE to
# /home/yeray/projects/atom-to-reactor.
CODE=/home/yeray/projects/atom-to-reactor-m1steps
OUT=/home/yeray/projects/atom-to-reactor/notebooks/results/fit8_M1-steps
PY=/home/yeray/miniforge3/envs/atom-to-reactor/bin/python
WORKERS=$(nproc)
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
        echo "Eighth fit, M1 tube by tube (setup: $OUT/fit_config.json), launched $(date '+%Y-%m-%d %H:%M:%S %Z') on $(hostname), $WORKERS workers"
        echo "code: $CODE, branch $(git rev-parse --abbrev-ref HEAD), commit $(git rev-parse HEAD) plus the uncommitted changes in code_vs_commit.patch (kinetics/, scripts/, tests/, data/)"
        echo "untracked code files: $(git ls-files --others --exclude-standard -- kinetics scripts tests | tr '\n' ' ')"
        echo "observables: $(md5sum "$DATA/data/lab_observables.json")"
        echo "snapshot (solvation): $(md5sum "$DATA/data/tank_api_snapshot_2026-10-08.json")"
        echo "snapshot (gas phase, geometries, shieldings; identical in the two files): $(md5sum "$DATA/data/tank_api_snapshot_2026-09-24.json")"
        echo "stages: registered, steps ${QUICK[*]}"
    } > "$OUT/provenance.txt"
else
    echo "[$(date +%T)] chain restarted at commit $(git rev-parse HEAD) (argument: ${1:-none}); finished jobs are read from their files" >> "$OUT/provenance.txt"
fi

run() {
    "$PY" scripts/run_fit.py "$@" "${QUICK[@]}" --workers "$WORKERS" --results "$OUT" >> "$OUT/fit8.out" 2>&1
    code=$?
    if [ $code -ne 0 ]; then
        echo "[$(date +%T)] CHAIN STOPPED: run_fit.py $* failed (exit $code), see $OUT/fit8.out" | tee -a "$OUT/fit8.out" "$OUT/run.log"
        exit 1
    fi
}

run registered
run steps
echo "[$(date +%T)] CHAIN FINISHED" | tee -a "$OUT/fit8.out" "$OUT/run.log"
