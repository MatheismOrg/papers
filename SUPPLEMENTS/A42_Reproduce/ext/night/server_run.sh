#!/bin/bash
# SA74 night, server side (Gary 23:24: "u can use my more powerful server also"). Runs every chain in jobs_server.txt
# (TAG L K THERM MEAS SEED A) to completion at low priority (nice 19), on all cores but two, in segments of SEG CPU seconds
# that save TAG.cfg, so progress survives a stop. Engine: regge_x86, a static build of regge_r2_big.c (no compiler needed).
# Stop: pkill -f sa74_night
cd "$(dirname "$0")"
export RG="$PWD/regge_x86"; chmod +x "$RG"
N=$(nproc); P=${P:-$(( N > 3 ? N - 2 : 1 ))}; SEG=${SEG:-600}
echo "$(date +%H:%M) start: $N cores, load $(cut -d' ' -f1-3 /proc/loadavg), running $P chains at a time, segments of $SEG s"
one() {
  tag=$1; L=$2; k=$3; therm=$4; meas=$5; seed=$6; a=$7
  while :; do
    th=0; me=0
    if [ -f $tag.dat ]; then th=$(( $(awk '$2==0' $tag.dat | wc -l) )); me=$(( $(awk '$2==1' $tag.dat | wc -l) )); fi
    if [ $me -ge $meas ]; then echo "$(date +%H:%M) $tag complete th=$th me=$me"; return; fi
    part=$(( $(ls $tag.part*.dat 2>/dev/null | wc -l) )); load=-; [ -f $tag.cfg ] && load=$tag.cfg
    nt=$(( therm > th ? therm - th : 0 )); nm=$(( meas - me ))
    MAXSEC=$SEG R2A=$a nice -n 19 "$RG" run $L $k 1 $nt $nm $((seed + 1000 * part)) $tag.part$part.dat 0.5 $load $tag.cfg 10 0 0 > /dev/null || { echo "$tag failed"; return; }
    grep -v '^#' $tag.part$part.dat >> $tag.dat; grep '^# done' $tag.part$part.dat | sed "s/^/# part$part /" >> $tag.log
    echo "$(date +%H:%M) $tag p$part th=$(( $(awk '$2==0' $tag.dat | wc -l) )) me=$(( $(awk '$2==1' $tag.dat | wc -l) ))"
  done
}
export -f one; export SEG
grep -v '^#' jobs_server.txt | xargs -P $P -L 1 bash -c 'one "$@"' _
echo "$(date +%H:%M) ALL DONE"
