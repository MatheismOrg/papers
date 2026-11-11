#!/bin/bash
# SA74, server side: runs the one job in jobs_s3.txt (TAG L K THERM MEAS SEED A) with the static multi-core engine regge_par_x86
# (PAR=$NTH threads, Wilson lines) in segments of SEG wall-clock seconds that save TAG.cfg. Started at idle priority (chrt -i 0) by
# start_server3.sh, so the other project and any other job keep their CPU. Stop: pkill -f "[s]a74_s3"
cd "$(dirname "$0")"
RG="$PWD/regge_par_x86"; SEG=${SEG:-600}; NTH=${NTH:-6}
read -r tag L k therm meas seed a < <(grep -v '^#' jobs_s3.txt | head -1)
echo "$(date +%H:%M) start $tag on $(nproc) cores, load $(cut -d' ' -f1-3 /proc/loadavg), $NTH threads, idle priority"
while :; do
  th=0; me=0
  if [ -f $tag.dat ]; then th=$(awk '$2==0' $tag.dat | wc -l); me=$(awk '$2==1' $tag.dat | wc -l); fi
  if [ $me -ge $meas ]; then echo "$(date +%H:%M) $tag complete th=$th me=$me"; break; fi
  part=$(ls $tag.part*.dat 2>/dev/null | wc -l); load=-; [ -f $tag.cfg ] && load=$tag.cfg
  nt=$(( therm > th ? therm - th : 0 )); nm=$(( meas - me ))
  PAR=$NTH LINES=1 MAXSEC=$SEG R2A=$a "$RG" run $L $k 1 $nt $nm $((seed + 1000 * part)) $tag.part$part.dat 0.5 $load $tag.cfg 10 0 0 > /dev/null || { echo "$(date +%H:%M) $tag failed"; break; }
  grep -v '^#' $tag.part$part.dat >> $tag.dat; grep '^# done' $tag.part$part.dat | sed "s/^/# part$part /" >> $tag.log
  echo "$(date +%H:%M) $tag p$part th=$(awk '$2==0' $tag.dat | wc -l) me=$(awk '$2==1' $tag.dat | wc -l)"
done
