#!/bin/bash
# SA74 runner 2 (Gary 03:24: "u can continue past 4am without limit"): runs the jobs of queue2.txt (TAG L K THERM MEAS SEED A [FROM])
# with the sequential-order engine regge_mac_seq (SEQ=1) and Wilson-line output (LINES=1), with no end time. It shares the cores with
# runner 1 (night_queue.sh): it starts a job only while fewer than (cores - 1) engines run in total, and before 05:30 only when runner 1
# has no startable job left. A file STOP2 in this folder makes every runner-2 job stop after its current segment (progress kept, claim
# released) and the runner exit. Stop everything at once: pkill -f sa74_night_regge
cd "$(dirname "$0")"
export RG2="$PWD/regge_mac_seq" SEG=${SEG:-600}
MAXP=${MAXP:-$(( $(sysctl -n hw.ncpu) - 1 ))}
count() { pgrep -f "$PWD/regge_mac" | wc -l | tr -d ' '; }          # engines of both runners
pending1() {                                                          # startable jobs left in runner 1's queue
  local n=0; [ -f queue.txt ] || { echo 0; return; }
  while read -r tag L k th me seed a from; do case "$tag" in ''|\#*) continue;; esac
    [ -d $tag.claim ] && continue; [ -n "$from" ] && [ ! -f $from.done ] && continue; n=$((n + 1)); done < queue.txt; echo $n; }
defer() { pgrep -f "[n]ight_queue\.sh" > /dev/null && [ "$(date +%H%M)" -lt 0530 ] && [ "$(pending1)" -gt 0 ]; }
run_job() {
  tag=$1; L=$2; k=$3; therm=$4; meas=$5; seed=$6; a=$7; from=$8
  [ -n "$from" ] && [ ! -f $tag.cfg ] && cp $from.cfg $tag.cfg
  while :; do
    if [ -f STOP2 ]; then rmdir $tag.claim 2>/dev/null; echo "$(date +%H:%M) $tag paused (STOP2)"; return; fi
    th=0; me=0
    if [ -f $tag.dat ]; then th=$(( $(awk '$2==0' $tag.dat | wc -l) )); me=$(( $(awk '$2==1' $tag.dat | wc -l) )); fi
    if [ $me -ge $meas ]; then touch $tag.done; echo "$(date +%H:%M) $tag complete th=$th me=$me"; return; fi
    part=$(( $(ls $tag.part*.dat 2>/dev/null | wc -l) )); load=-; [ -f $tag.cfg ] && load=$tag.cfg
    nt=$(( therm > th ? therm - th : 0 )); nm=$(( meas - me ))
    SEQ=1 LINES=1 MAXSEC=$SEG R2A=$a "$RG2" run $L $k 1 $nt $nm $((seed + 1000 * part)) $tag.part$part.dat 0.5 $load $tag.cfg 10 0 0 > /dev/null || { echo "$(date +%H:%M) $tag failed"; touch $tag.failed; return; }
    grep -v '^#' $tag.part$part.dat >> $tag.dat; grep '^# done' $tag.part$part.dat | sed "s/^/# part$part /" >> $tag.log
    echo "$(date +%H:%M) $tag p$part th=$(( $(awk '$2==0' $tag.dat | wc -l) )) me=$(( $(awk '$2==1' $tag.dat | wc -l) ))"
  done
}
echo "$(date +%H:%M) runner 2: up to $MAXP engine processes in total, segments of $SEG s, sequential engine with Wilson lines"
while [ ! -f STOP2 ]; do
  r1=$(count); sleep 3; r2=$(count); running=$(( r1 > r2 ? r1 : r2 )); started=0
  if [ $running -lt $MAXP ] && ! defer; then
    while read -r tag L k therm meas seed a from; do
      case "$tag" in ''|\#*) continue;; esac
      [ -d $tag.claim ] && continue; [ -f $tag.done ] && continue
      [ -n "$from" ] && [ ! -f $from.done ] && continue
      mkdir $tag.claim 2>/dev/null || continue
      echo "$(date +%H:%M) start $tag (engines running before: $running)"
      run_job $tag $L $k $therm $meas $seed $a $from < /dev/null &
      started=1; break
    done < queue2.txt
  fi
  if [ $started = 1 ]; then sleep 5; else sleep 30; fi
done
echo "$(date +%H:%M) runner 2: STOP2 found; jobs stop after their current segments"; wait
