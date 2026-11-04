#!/usr/bin/env bash
# Third fit (docs/plan-fit_improvement.md): the stages one after another; a failed stage stops the chain.
# Skipped on purpose: 'refit' (repaired a solver setting of the first run) and 'throughput'.
cd /home/yeray/projects/atom-to-reactor || exit 1
PY=/home/yeray/miniforge3/envs/atom-to-reactor/bin/python
OUT=notebooks/results/05_fit3
for stage in registered recovery baselines extended night traces predictions; do
    "$PY" scripts/run_fit.py "$stage" --workers 12 --results "$OUT" >> "$OUT/fit3.out" 2>&1
    code=$?
    if [ $code -ne 0 ]; then
        echo "[$(date +%T)] CHAIN STOPPED: stage $stage failed (exit $code), see $OUT/fit3.out" | tee -a "$OUT/fit3.out" "$OUT/run.log"
        exit 1
    fi
done
echo "[$(date +%T)] CHAIN FINISHED: all stages done" | tee -a "$OUT/fit3.out" "$OUT/run.log"
