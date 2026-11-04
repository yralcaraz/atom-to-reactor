#!/usr/bin/env bash
# The eighth fit on a Google Cloud VM (2026-10-23): n4-highcpu-64 Spot (64 vCPU, the project's quota over all
# regions; 128 GB) in europe-west1-c, project atom-to-reactor-512413. The VM gets the same paths as the laptop
# (conda env, worktree, main .git for the worktree, 05_fit8, the fits of 05_fit7 it starts from, lab data), so
# run_chain.sh runs there unchanged and takes every core (nproc).
#
#   bash cloud.sh create    creates the VM (Spot; Google stops it after 12 h at most, the disk is kept)
#   bash cloud.sh setup     installs tmux/rsync/git and copies env, code, data, 05_fit8 and the fits of 05_fit7
#   bash cloud.sh smoke     the whole chain with --quick in 05_fit8_smoke, waits for it and shows its end
#   bash cloud.sh start     the run (tmux session fit8 on the VM)
#   bash cloud.sh status    VM state, the steps done and the last lines of run.log on the VM
#   bash cloud.sh resume    after Google stopped the VM (Spot preemption or 12 h): starts it and the chain again
#   bash cloud.sh back      copies 05_fit8 (and 05_fit8_smoke) from the VM to the laptop
#   bash cloud.sh delete    deletes the VM and its disk (after back)
#   bash cloud.sh budget    budget alert of 10 (billing account currency) on the project: e-mail only, stops nothing
#
# Data on the VM: lab_observables.json and the two tank snapshots (data/ of the private folder), nothing of the raw NMR.
set -euo pipefail
PROJECT=atom-to-reactor-512413
ZONE=europe-west1-c
VM=fit8
HOST=$VM.$ZONE.$PROJECT                     # alias written by `gcloud compute config-ssh`
BILLING=01B250-532505-3FE634
ENV=/home/yeray/miniforge3/envs/atom-to-reactor
REPO=/home/yeray/projects/atom-to-reactor
CODE=/home/yeray/projects/atom-to-reactor-fit8
OUT=$REPO/notebooks/results/05_fit8
FROM=$REPO/notebooks/results/05_fit7/fits
DATA=/home/yeray/data/atom_to_reactor_private/data
G="gcloud --project=$PROJECT"

ssh_vm() { ssh -o ServerAliveInterval=60 "$HOST" "$@"; }

case "${1:-}" in
create)
    $G compute instances create $VM --zone=$ZONE --machine-type=n4-highcpu-64 \
        --provisioning-model=SPOT --instance-termination-action=STOP --max-run-duration=12h \
        --image-family=ubuntu-2404-lts-amd64 --image-project=ubuntu-os-cloud \
        --boot-disk-size=20GB --boot-disk-type=hyperdisk-balanced \
        --no-service-account --no-scopes
    sleep 30
    $G compute ssh $VM --zone=$ZONE --quiet --command=true    # creates the SSH key on the first use
    $G compute config-ssh
    ;;
setup)
    ssh_vm "sudo apt-get -qq update && sudo apt-get -qq install -y rsync tmux git
        mkdir -p $(dirname $ENV) $FROM $DATA"
    rsync -az --info=progress2 $ENV/ "$HOST:$ENV/"
    rsync -az $REPO/.git/ "$HOST:$REPO/.git/"
    rsync -az --exclude=__pycache__ $CODE/ "$HOST:$CODE/"
    rsync -az $DATA/lab_observables.json $DATA/tank_api_snapshot_2026-09-24.json \
        $DATA/tank_api_snapshot_2026-10-08.json "$HOST:$DATA/"
    rsync -az $FROM/ "$HOST:$FROM/"
    rsync -az $OUT/ "$HOST:$OUT/"
    ssh_vm "nproc; free -g | head -2; $ENV/bin/python -c 'import numpy, scipy, pandas; print(\"env ok\")'"
    ;;
smoke)
    ssh_vm "cd $REPO && time bash $OUT/run_chain.sh smoke; tail -n 12 ${OUT}_smoke/run.log | cut -c1-160"
    ;;
start)
    if ssh_vm "tmux has-session -t fit8 2>/dev/null"; then echo "the chain already runs on the VM"; exit 1; fi
    ssh_vm "tmux new -d -s fit8 bash $OUT/run_chain.sh"
    echo "[$(date +%T)] chain started on the VM; follow it with: bash cloud.sh status"
    ;;
status)
    $G compute instances describe $VM --zone=$ZONE --format="value(status)"
    ssh_vm "uptime; grep -h 'STEP DONE\|CHAIN' $OUT/run.log; tail -n 5 $OUT/run.log | cut -c1-160"
    ;;
resume)
    $G compute instances start $VM --zone=$ZONE
    $G compute config-ssh                             # the external IP changes at every start
    sleep 30
    ssh_vm "tmux new -d -s fit8 bash $OUT/run_chain.sh"
    ;;
back)
    rsync -az "$HOST:$OUT/" $OUT/
    rsync -az "$HOST:${OUT}_smoke/" ${OUT}_smoke/ || true
    tail -n 3 $OUT/run.log
    ;;
delete)
    $G compute instances delete $VM --zone=$ZONE --delete-disks=all
    ;;
budget)
    gcloud billing budgets create --billing-account=$BILLING --display-name="fit8 VM" --budget-amount=10 \
        --filter-projects=projects/$PROJECT --threshold-rule=percent=0.5 --threshold-rule=percent=1.0
    ;;
*)
    sed -n '2,17p' "$0"; exit 1
    ;;
esac
