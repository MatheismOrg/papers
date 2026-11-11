#!/bin/bash
# SA74 night: native macOS runner (Gary 23:11: "does it help u to use my local gpu or cpu? if so u can").
# Runs every chain in jobs_native.txt (TAG L K THERM MEAS SEED A) to completion, P chains at a time, in segments of SEG CPU
# seconds; each segment loads and saves TAG.cfg and appends its rows to TAG.dat, so progress survives a stop.
# Start (paste into Terminal): bash ~/repos/gr-qm-tests/sa74_night_regge/start_native.sh
# Stop: pkill -f sa74_night_regge
cd "$(dirname "$0")"
while grep -q busy vm.lock 2>/dev/null; do sleep 5; done   # let a running in-app chunk finish first
clang -O2 -o regge_mac regge_r2_big.c || exit 1
export RG="$PWD/regge_mac"
P=${P:-$(( $(sysctl -n hw.ncpu) - 1 ))}
SEG=${SEG:-600}
echo "$(date +%H:%M) start: $(sysctl -n hw.ncpu) cores, running $P chains at a time, segments of $SEG s"
one() {
  tag=$1; L=$2; k=$3; therm=$4; meas=$5; seed=$6; a=$7
  while :; do
    th=0; me=0
    if [ -f $tag.dat ]; then th=$(( $(awk '$2==0' $tag.dat | wc -l) )); me=$(( $(awk '$2==1' $tag.dat | wc -l) )); fi
    if [ $me -ge $meas ]; then echo "$(date +%H:%M) $tag complete th=$th me=$me"; return; fi
    part=$(( $(ls $tag.part*.dat 2>/dev/null | wc -l) )); load=-; [ -f $tag.cfg ] && load=$tag.cfg
    nt=$(( therm > th ? therm - th : 0 )); nm=$(( meas - me ))
    MAXSEC=$SEG R2A=$a "$RG" run $L $k 1 $nt $nm $((seed + 1000 * part)) $tag.part$part.dat 0.5 $load $tag.cfg 10 0 0 > /dev/null || { echo "$tag failed"; return; }
    grep -v '^#' $tag.part$part.dat >> $tag.dat; grep '^# done' $tag.part$part.dat | sed "s/^/# part$part /" >> $tag.log
    echo "$(date +%H:%M) $tag p$part th=$(( $(awk '$2==0' $tag.dat | wc -l) )) me=$(( $(awk '$2==1' $tag.dat | wc -l) ))"
  done
}
export -f one; export SEG
grep -v '^#' jobs_native.txt | xargs -P $P -L 1 bash -c 'one "$@"' _
echo "$(date +%H:%M) ALL DONE"
