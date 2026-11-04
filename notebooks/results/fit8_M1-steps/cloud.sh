#!/usr/bin/env bash
# The eighth fit, M1 tube by tube, on a Google Cloud VM (2026-10-23): n4-highcpu-16 Spot (16 vCPU, 32 GB),
# europe-west1-c, project atom-to-reactor-512413. Spot: if Google takes the VM back it is deleted, not stopped (the
# run takes minutes: create it again and start over). The VM gets
# the same paths as the laptop (conda env, worktree, main .git for the worktree, fit8_M1-steps, lab data), so
# run_chain.sh runs there unchanged and takes every core (nproc). Adapted from 05_fit8/cloud.sh.
#
#   bash cloud.sh create    creates the VM; Google deletes it with its disk after 3 h, or when it takes a Spot VM back
#   bash cloud.sh setup     installs tmux/rsync/git and copies env, code, data and fit8_M1-steps
#   bash cloud.sh smoke     the chain with --quick in fit8_M1-steps_smoke, waits for it and shows its end
#   bash cloud.sh start     the run (tmux session fit8 on the VM)
#   bash cloud.sh status    VM state and the last lines of run.log on the VM
#   bash cloud.sh back      copies fit8_M1-steps (and its smoke) from the VM to the laptop
#   bash cloud.sh delete    deletes the VM and its disk (after back), then checks that no VM and no disk is left
#
# Data on the VM: lab_observables.json and the two tank snapshots (data/ of the private folder), nothing of the raw NMR.
# They live on the VM's boot disk only, which is deleted with the VM.
set -euo pipefail
PROJECT=atom-to-reactor-512413
ZONE=europe-west1-c
VM=fit8-m1steps
HOST=$VM.$ZONE.$PROJECT                     # alias written by `gcloud compute config-ssh`
ENV=/home/yeray/miniforge3/envs/atom-to-reactor
REPO=/home/yeray/projects/atom-to-reactor
CODE=/home/yeray/projects/atom-to-reactor-m1steps
OUT=$REPO/notebooks/results/fit8_M1-steps
DATA=/home/yeray/data/atom_to_reactor_private/data
G="gcloud --project=$PROJECT"

ssh_vm() { ssh -o ServerAliveInterval=60 "$HOST" "$@"; }

case "${1:-}" in
create)
    $G compute instances create $VM --zone=$ZONE --machine-type=n4-highcpu-16 \
        --provisioning-model=SPOT --max-run-duration=3h --instance-termination-action=DELETE \
        --image-family=ubuntu-2404-lts-amd64 --image-project=ubuntu-os-cloud \
        --boot-disk-size=20GB --boot-disk-type=hyperdisk-balanced --boot-disk-auto-delete \
        --no-service-account --no-scopes
    sleep 30
    $G compute ssh $VM --zone=$ZONE --quiet --command=true    # creates the SSH key on the first use
    $G compute config-ssh
    ;;
setup)
    ssh_vm "sudo apt-get -qq update && sudo apt-get -qq install -y rsync tmux git
        mkdir -p $(dirname $ENV) $DATA $REPO $OUT"
    rsync -az --info=progress2 $ENV/ "$HOST:$ENV/"
    rsync -az $REPO/.git/ "$HOST:$REPO/.git/"
    rsync -az --exclude=__pycache__ $CODE/ "$HOST:$CODE/"
    rsync -az $DATA/lab_observables.json $DATA/tank_api_snapshot_2026-09-24.json \
        $DATA/tank_api_snapshot_2026-10-08.json "$HOST:$DATA/"
    rsync -az $OUT/fit_config.json $OUT/run_chain.sh "$HOST:$OUT/"
    ssh_vm "nproc; free -g | head -2; $ENV/bin/python -c 'import numpy, scipy, pandas; print(\"env ok\")'"
    ;;
smoke)
    ssh_vm "cd $REPO && time bash $OUT/run_chain.sh smoke; tail -n 8 ${OUT}_smoke/run.log | cut -c1-160"
    ;;
start)
    if ssh_vm "tmux has-session -t fit8 2>/dev/null"; then echo "the chain already runs on the VM"; exit 1; fi
    ssh_vm "tmux new -d -s fit8 bash $OUT/run_chain.sh"
    echo "[$(date +%T)] chain started on the VM; follow it with: bash cloud.sh status"
    ;;
status)
    $G compute instances describe $VM --zone=$ZONE --format="value(status)"
    ssh_vm "uptime; grep -h 'CHAIN' $OUT/run.log || true; tail -n 5 $OUT/run.log | cut -c1-160"
    ;;
back)
    rsync -az "$HOST:$OUT/" $OUT/
    rsync -az "$HOST:${OUT}_smoke/" ${OUT}_smoke_cloud/ || true
    tail -n 3 $OUT/run.log
    ;;
delete)
    $G compute instances delete $VM --zone=$ZONE --delete-disks=all --quiet || true
    echo "VMs left in the project:";   $G compute instances list --format="value(name,zone,status)"
    echo "disks left in the project:"; $G compute disks list --format="value(name,zone)"
    ;;
*)
    sed -n '2,16p' "$0"; exit 1
    ;;
esac
