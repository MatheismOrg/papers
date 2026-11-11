#!/bin/bash
# SA74 runner 3 (plan A5 and later): runs the jobs of queue3.txt (TAG L K THERM MEAS SEED A [FROM]) one at a time with the multi-core
# engine regge_mac_par (PAR=$NTH threads, Wilson lines), in segments of SEG wall-clock seconds that save TAG.cfg, so progress survives a
# stop. It starts a job only while no other regge_mac engine runs, so it waits for runner 2's chains. A file STOP3 in this folder stops
# it after the current segment. Stop everything at once: pkill -f sa74_night_regge
cd "$(dirname "$0")"
RG3="$PWD/regge_mac_par"; SEG=${SEG:-600}; NTH=${NTH:-8}
busy() { pgrep -f "$PWD/regge_mac" > /dev/null; }
echo "$(date +%H:%M) runner 3: one job at a time with $NTH threads; it waits while another engine runs"
while [ ! -f STOP3 ]; do
  job=""
  while read -r tag L k therm meas seed a from; do
    case "$tag" in ''|\#*) continue;; esac
    [ -f $tag.done ] && continue; [ -f $tag.failed ] && continue
    [ -n "$from" ] && [ ! -f $from.done ] && continue
    job="$tag $L $k $therm $meas $seed $a $from"; break
  done < queue3.txt
  if [ -z "$job" ] || busy; then sleep 60; continue; fi
  set -- $job; tag=$1; L=$2; k=$3; therm=$4; meas=$5; seed=$6; a=$7; from=$8
  [ -n "$from" ] && [ ! -f $tag.cfg ] && cp $from.cfg $tag.cfg
  echo "$(date +%H:%M) start $tag"
  while [ ! -f STOP3 ]; do
    th=0; me=0
    if [ -f $tag.dat ]; then th=$(awk '$2==0' $tag.dat | wc -l | tr -d ' '); me=$(awk '$2==1' $tag.dat | wc -l | tr -d ' '); fi
    if [ $me -ge $meas ]; then touch $tag.done; echo "$(date +%H:%M) $tag complete th=$th me=$me"; break; fi
    part=$(ls $tag.part*.dat 2>/dev/null | wc -l | tr -d ' '); load=-; [ -f $tag.cfg ] && load=$tag.cfg
    nt=$(( therm > th ? therm - th : 0 )); nm=$(( meas - me ))
    PAR=$NTH LINES=1 MAXSEC=$SEG R2A=$a "$RG3" run $L $k 1 $nt $nm $((seed + 1000 * part)) $tag.part$part.dat 0.5 $load $tag.cfg 10 0 0 > /dev/null || { echo "$(date +%H:%M) $tag failed"; touch $tag.failed; break; }
    grep -v '^#' $tag.part$part.dat >> $tag.dat; grep '^# done' $tag.part$part.dat | sed "s/^/# part$part /" >> $tag.log
    echo "$(date +%H:%M) $tag p$part th=$(awk '$2==0' $tag.dat | wc -l | tr -d ' ') me=$(awk '$2==1' $tag.dat | wc -l | tr -d ' ')"
  done
done
echo "$(date +%H:%M) runner 3: STOP3 found; stopped"
