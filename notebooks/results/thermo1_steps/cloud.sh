#!/usr/bin/env bash
# The thermodynamic fit thermo1 on a Google Cloud VM (prepared 2026-10-29): n4-highcpu-64 Spot (64 vCPU, 128 GB),
# europe-west1-c, project atom-to-reactor-512413. Adapted from fit8_M1-steps/cloud.sh. Spot: if Google takes the VM back
# it is deleted with its disk, not stopped (create it again and start over). The VM gets the same paths as the laptop
# (conda env, repository working tree, lab data), so run_parallel.sh runs there unchanged: 6 parts (ξ × chain) × 10 workers.
#
#   bash cloud.sh create    creates the VM; Google deletes it with its disk after 3 h, or when it takes a Spot VM back
#   bash cloud.sh setup     installs tmux/rsync, copies env, code (working tree, uncommitted fit code included), data and
#                           the three input files of thermo1_steps; runs tests/test_thermo_fit.py there
#   bash cloud.sh smoke     the parallel run with --quick in thermo1_steps_smoke, waits for it and shows its end
#   bash cloud.sh start     the run (tmux session thermo1 on the VM)
#   bash cloud.sh status    VM state, load and the last lines of run.log on the VM
#   bash cloud.sh back      copies thermo1_steps (and its smoke, as thermo1_steps_smoke_cloud) from the VM to the laptop,
#                           then compares both by checksum: "BACK VERIFIED" only if every file on the VM is here unchanged
#   bash cloud.sh delete    wipes the lab data and results on the VM, deletes the VM and its disk, removes the SSH alias
#   bash cloud.sh verify    after delete: the results complete on the laptop, and nothing left in the project (VMs,
#                           disks, snapshots, images, machine images, reserved IPs, buckets): "ALL CLEAN" or what is left
#
# Data on the VM: lab_observables.json, lab_sample_folders.json, the two tank snapshots, and the shares of the spectra at
# equilibrium (spectra_shares.csv, equilibrium_points.csv). Nothing of the raw NMR: the stage 'points' is never run there.
# They live on the VM's boot disk only, which is deleted with the VM.
set -euo pipefail
PROJECT=atom-to-reactor-512413
ZONE=europe-west1-c
VM=thermo1
MACHINE=n4-highcpu-64
HOST=$VM.$ZONE.$PROJECT                     # alias written by `gcloud compute config-ssh`
ENV=/home/yeray/miniforge3/envs/atom-to-reactor
REPO=/home/yeray/projects/atom-to-reactor
OUT=$REPO/notebooks/results/thermo1_steps
SMOKE=${OUT}_smoke
DATA=/home/yeray/data/atom_to_reactor_private/data
PY=$ENV/bin/python
G="gcloud --project=$PROJECT"

ssh_vm() { ssh -o ServerAliveInterval=60 "$HOST" "$@"; }

case "${1:-}" in
create)
    $G compute instances create $VM --zone=$ZONE --machine-type=$MACHINE \
        --provisioning-model=SPOT --max-run-duration=3h --instance-termination-action=DELETE \
        --image-family=ubuntu-2404-lts-amd64 --image-project=ubuntu-os-cloud \
        --boot-disk-size=20GB --boot-disk-type=hyperdisk-balanced --boot-disk-auto-delete \
        --no-service-account --no-scopes
    sleep 30
    $G compute ssh $VM --zone=$ZONE --quiet --command=true    # creates the SSH key on the first use
    $G compute config-ssh
    ;;
setup)
    ssh_vm "sudo apt-get -qq update && sudo apt-get -qq install -y rsync tmux
        mkdir -p $(dirname $ENV) $DATA $OUT"
    rsync -az --info=progress2 $ENV/ "$HOST:$ENV/"
    rsync -az --exclude=.git --exclude=__pycache__ --exclude=notebooks/results $REPO/ "$HOST:$REPO/"
    rsync -az $DATA/lab_observables.json $DATA/lab_sample_folders.json $DATA/tank_api_snapshot_2026-09-24.json \
        $DATA/tank_api_snapshot_2026-10-08.json "$HOST:$DATA/"
    rsync -az $OUT/fit_config.json $OUT/spectra_shares.csv $OUT/equilibrium_points.csv $OUT/run_parallel.sh "$HOST:$OUT/"
    ssh_vm "nproc; free -g | head -2; cd $REPO && $PY tests/test_thermo_fit.py"
    ;;
smoke)
    ssh_vm "rm -rf $SMOKE && mkdir -p $SMOKE && cp $OUT/fit_config.json $OUT/spectra_shares.csv $OUT/equilibrium_points.csv $SMOKE/
        cd $REPO && time bash $OUT/run_parallel.sh $SMOKE --quick; tail -n 3 $SMOKE/run.log | cut -c1-200"
    ;;
start)
    if ssh_vm "tmux has-session -t thermo1 2>/dev/null"; then echo "the run already runs on the VM"; exit 1; fi
    if ssh_vm "test -d $OUT/steps"; then echo "$OUT/steps exists on the VM: the run was started before (it resumes from there)"; fi
    ssh_vm "tmux new -d -s thermo1 bash $OUT/run_parallel.sh $OUT"
    echo "[$(date +%T)] run started on the VM; follow it with: bash cloud.sh status"
    ;;
status)
    $G compute instances describe $VM --zone=$ZONE --format="value(status)"
    ssh_vm "uptime; ls $OUT/steps 2>/dev/null | wc -l | sed 's/^/fits on disk (of 36): /'; grep -h 'PARALLEL RUN DONE' $OUT/run.log || true
        tail -n 5 $OUT/run.log | cut -c1-200"
    ;;
back)
    rsync -az "$HOST:$OUT/" $OUT/
    rsync -az "$HOST:$SMOKE/" ${OUT}_smoke_cloud/ || true
    # a second pass by checksum that would copy nothing: every file of the VM is on the laptop, byte for byte
    left=$(rsync -anc --itemize-changes "$HOST:$OUT/" $OUT/ | grep -v '^\.d' || true)
    remote=$(ssh_vm "find $OUT -type f | wc -l"); here=$(find $OUT -type f | wc -l)
    echo "files: $remote on the VM, $here here; fits here: $(find $OUT/steps -name "*.json" 2>/dev/null | wc -l) of 36"
    if [ -z "$left" ] && [ "$here" -ge "$remote" ]; then echo "BACK VERIFIED"; else echo "NOT VERIFIED, differs:"; echo "$left"; exit 1; fi
    tail -n 3 $OUT/run.log | cut -c1-200
    ;;
delete)
    # the lab data and results off the VM first (the disk goes with the VM anyway; this is a second guard)
    ssh_vm "rm -rf $DATA $OUT $SMOKE && echo 'lab data and results wiped on the VM'" || echo "VM not reachable: nothing wiped by hand"
    $G compute instances delete $VM --zone=$ZONE --delete-disks=all --quiet || true
    $G compute config-ssh --remove >/dev/null || true          # the SSH alias of the VM out of ~/.ssh/config
    echo "VMs left in the project:";   $G compute instances list --format="value(name,zone,status)"
    echo "disks left in the project:"; $G compute disks list --format="value(name,zone)"
    ;;
verify)
    set +e                                                  # every check runs, a failing one only marks NOT CLEAN
    echo "== results on the laptop ($OUT)"
    ok=1
    n=$(find $OUT/steps -name "*.json" 2>/dev/null | wc -l); echo "fits: $n of 36"; [ "$n" -eq 36 ] || ok=0
    for f in steps_summary.csv steps_scores.csv steps_intervals.csv steps_energies.csv replicates.csv; do
        if [ -s $OUT/$f ]; then echo "ok      $f"; else echo "MISSING $f"; ok=0; fi
    done
    done_line=$(grep -h 'PARALLEL RUN DONE' $OUT/run.log | tail -n 1); echo "${done_line:-no PARALLEL RUN DONE line}"
    echo "$done_line" | grep -q 'tables=0' && ! echo "$done_line" | grep -qE '=[1-9]' || ok=0
    echo "== Google Cloud project $PROJECT (each list must be empty)"
    for what in "compute instances" "compute disks" "compute snapshots" "compute machine-images" "compute addresses"; do
        out=$($G $what list --format="value(name)" 2>&1); echo "$what: ${out:-none}"; [ -z "$out" ] || ok=0
    done
    out=$($G compute images list --no-standard-images --format="value(name)" 2>&1); echo "compute images (own): ${out:-none}"; [ -z "$out" ] || ok=0
    out=$($G storage buckets list --format="value(name)" 2>&1); echo "storage buckets: ${out:-none}"; [ -z "$out" ] || ok=0
    if grep -q "$HOST" ~/.ssh/config 2>/dev/null; then echo "SSH alias still in ~/.ssh/config"; ok=0; else echo "SSH alias: removed"; fi
    if [ $ok -eq 1 ]; then echo "ALL CLEAN: results complete here, nothing left in the cloud"; else echo "NOT CLEAN: see above"; exit 1; fi
    ;;
*)
    sed -n '2,20p' "$0"; exit 1
    ;;
esac
