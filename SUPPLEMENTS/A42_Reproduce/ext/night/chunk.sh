#!/bin/bash
# One time-boxed chunk (MAXSEC CPU seconds, default 140) for each chain in jobs.txt (TAG L K THERM MEAS SEED A), all in parallel.
# Each chunk loads TAG.cfg if present, runs the remaining thermalization (phase 0) and measurement (phase 1) sweeps until the
# time box ends after a whole sweep, saves TAG.cfg and appends its rows to TAG.dat; its '# done' line goes to TAG.log.
cd "$(dirname "$0")"
[ -f night_native.log ] && { echo "native runner active; in-app chunks stop"; exit 0; }
echo busy > vm.lock; trap 'echo free > vm.lock' EXIT TERM INT
[ -x regge_r2_big ] || gcc -O2 -o regge_r2_big regge_r2_big.c -lm
MS=${MAXSEC:-140}
run_one() {
  tag=$1; L=$2; k=$3; therm=$4; meas=$5; seed=$6; a=$7
  th_done=0; me_done=0
  if [ -f $tag.dat ]; then th_done=$(awk '$2==0' $tag.dat | wc -l); me_done=$(awk '$2==1' $tag.dat | wc -l); fi
  [ "$me_done" -ge "$meas" ] && { echo "$tag complete th=$th_done me=$me_done"; return; }
  part=$(ls $tag.part*.dat 2>/dev/null | wc -l); load=-; [ -f $tag.cfg ] && load=$tag.cfg
  nt=$(( therm > th_done ? therm - th_done : 0 )); nm=$(( meas - me_done ))
  MAXSEC=$MS R2A=$a ./regge_r2_big run $L $k 1 $nt $nm $((seed + 1000 * part)) $tag.part$part.dat 0.5 $load $tag.cfg 10 0 0 > /dev/null
  grep -v '^#' $tag.part$part.dat >> $tag.dat; grep '^# done' $tag.part$part.dat | sed "s/^/# part$part /" >> $tag.log
  echo "$tag p$part th=$(awk '$2==0' $tag.dat | wc -l) me=$(awk '$2==1' $tag.dat | wc -l) R_last=$(tail -1 $tag.dat | awk '{print $5}')"
}
mapfile -t lines < <(grep -v '^#' jobs.txt)
for l in "${lines[@]}"; do run_one $l & done; wait
