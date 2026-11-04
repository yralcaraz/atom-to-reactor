#!/usr/bin/env bash
# Sixth fit: waits for mains power, runs every stage on a tiny search in a scratch folder (1 day only, one recovery
# data set per case), and starts the chain (run_chain.sh) only if all of them pass. Nothing is written to this
# folder before the chain starts, except this script's log.
# Run inside tmux: tmux new -d -s fit6 bash notebooks/results/05_fit6/start_when_plugged.sh
cd /home/yeray/projects/atom-to-reactor || exit 1
PY=/home/yeray/miniforge3/envs/atom-to-reactor/bin/python
OUT=notebooks/results/05_fit6
SMOKE=/tmp/fit6_smoke
LOG=$OUT/start_when_plugged.log
echo "[$(date '+%F %T')] waiting for mains power (battery $(cat /sys/class/power_supply/BAT0/capacity) %)" >> $LOG
until [ "$(cat /sys/class/power_supply/AC/online)" = 1 ]; do sleep 30; done
rm -rf "$SMOKE" && mkdir -p "$SMOKE"
"$PY" - "$OUT/fit_config.json" "$SMOKE/fit_config.json" <<'PYEOF'
import json, sys
config = json.load(open(sys.argv[1]))
config['scenarios'] = ['middle']
for case in config['recovery']['cases'].values():
    case['seeds'] = case['seeds'][:1]
json.dump(config, open(sys.argv[2], 'w'), indent=1)
PYEOF
echo "[$(date '+%F %T')] on mains power; smoke run of every stage in $SMOKE (1 day, tiny search)" >> $LOG
for st in registered recovery baselines extended night traces predictions; do
    "$PY" scripts/run_fit.py $st --quick --workers 12 --results "$SMOKE" > "$SMOKE/out_$st.txt" 2>&1 || {
        echo "[$(date '+%F %T')] SMOKE FAILED at stage $st: chain NOT started (see $SMOKE/out_$st.txt)" >> $LOG; exit 1; }
done
echo "[$(date '+%F %T')] smoke passed; chain started" >> $LOG
exec bash $OUT/run_chain.sh
