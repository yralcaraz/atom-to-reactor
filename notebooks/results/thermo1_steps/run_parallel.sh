#!/usr/bin/env bash
# The thermodynamic fit thermo1 in parallel: one process per ξ (1, 2, 3) and chain (a: probe, b: glovebox), WORKERS each,
# then one run without --xi/--chain that finds every fit on disk and writes the tables (seconds).
#
#   bash run_parallel.sh <results dir> [--quick]       (from anywhere; WORKERS=10 by default: 6 × 10 on a 64-vCPU VM)
#
# Each part logs to <results dir>/part_xi<ξ>_<chain>.out and, like the whole run, to run.log; the last line written to
# run.log is "PARALLEL RUN DONE" with the exit codes.
set -uo pipefail
REPO=/home/yeray/projects/atom-to-reactor
PY=/home/yeray/miniforge3/envs/atom-to-reactor/bin/python
OUT=$(realpath "$1"); shift
WORKERS=${WORKERS:-10}
cd "$REPO"
pids=(); names=()
for xi in 1 2 3; do
    for chain in a b; do
        $PY scripts/run_thermo_fit.py steps --results "$OUT" --xi $xi --chain $chain --workers "$WORKERS" "$@" \
            > "$OUT/part_xi${xi}_${chain}.out" 2>&1 &
        pids+=($!); names+=("xi${xi}_${chain}")
    done
done
codes=()
for i in "${!pids[@]}"; do
    wait "${pids[$i]}"; codes+=("${names[$i]}=$?")
done
$PY scripts/run_thermo_fit.py steps --results "$OUT" --workers "$WORKERS" "$@" > "$OUT/part_tables.out" 2>&1
codes+=("tables=$?")
echo "[$(date +%T)] PARALLEL RUN DONE: ${codes[*]}" >> "$OUT/run.log"
