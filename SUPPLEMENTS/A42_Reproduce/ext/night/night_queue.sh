#!/bin/bash
# SA74 night: a refillable job queue for the Mac's cores (Gary 00:06: "u can also use my local mac cpus").
# Every cycle, while fewer than (cores - 1) engine processes run on this Mac (night_native.sh's included), it starts the next job
# in queue.txt (TAG L K THERM MEAS SEED A [FROM]) that is not yet claimed and whose FROM job, if any, is done; a job with FROM
# starts from FROM's final configuration. Jobs run in segments of SEG CPU seconds that save TAG.cfg; a finished job leaves
# TAG.done. Claude can append jobs to queue.txt at any time. Stops starting jobs at 05:30. Stop everything: pkill -f sa74_night_regge
cd "$(dirname "$0")"
[ -x regge_mac ] || clang -O2 -o regge_mac regge_r2_big.c || exit 1
export RG="$PWD/regge_mac" SEG=${SEG:-600}
MAXP=${MAXP:-$(( $(sysctl -n hw.ncpu) - 1 ))}
count() { pgrep -f "$RG" | wc -l | tr -d ' '; }
run_job() {
  tag=$1; L=$2; k=$3; therm=$4; meas=$5; seed=$6; a=$7; from=$8
  [ -n "$from" ] && [ ! -f $tag.cfg ] && cp $from.cfg $tag.cfg
  while :; do
    th=0; me=0
    if [ -f $tag.dat ]; then th=$(( $(awk '$2==0' $tag.dat | wc -l) )); me=$(( $(awk '$2==1' $tag.dat | wc -l) )); fi
    if [ $me -ge $meas ]; then touch $tag.done; echo "$(date +%H:%M) $tag complete th=$th me=$me"; return; fi
    part=$(( $(ls $tag.part*.dat 2>/dev/null | wc -l) )); load=-; [ -f $tag.cfg ] && load=$tag.cfg
    nt=$(( therm > th ? therm - th : 0 )); nm=$(( meas - me ))
    MAXSEC=$SEG R2A=$a "$RG" run $L $k 1 $nt $nm $((seed + 1000 * part)) $tag.part$part.dat 0.5 $load $tag.cfg 10 0 0 > /dev/null || { echo "$(date +%H:%M) $tag failed"; touch $tag.failed; return; }
    grep -v '^#' $tag.part$part.dat >> $tag.dat; grep '^# done' $tag.part$part.dat | sed "s/^/# part$part /" >> $tag.log
    echo "$(date +%H:%M) $tag p$part th=$(( $(awk '$2==0' $tag.dat | wc -l) )) me=$(( $(awk '$2==1' $tag.dat | wc -l) ))"
  done
}
echo "$(date +%H:%M) queue runner: up to $MAXP engine processes, segments of $SEG s"
while [ "$(date +%H%M)" -lt 0530 ] || [ "$(date +%H)" -ge 12 ]; do
  r1=$(count); sleep 3; r2=$(count); running=$(( r1 > r2 ? r1 : r2 )); started=0
  if [ $running -lt $MAXP ]; then
    while read -r tag L k therm meas seed a from; do
      case "$tag" in ''|\#*) continue;; esac
      [ -d $tag.claim ] && continue
      [ -n "$from" ] && [ ! -f $from.done ] && continue
      mkdir $tag.claim 2>/dev/null || continue
      echo "$(date +%H:%M) start $tag (engines running before: $running)"
      run_job $tag $L $k $therm $meas $seed $a $from < /dev/null &
      started=1; break
    done < queue.txt
  fi
  if [ $started = 1 ]; then sleep 5; else sleep 30; fi
done
echo "$(date +%H:%M) queue runner stops starting jobs; running jobs finish"; wait
