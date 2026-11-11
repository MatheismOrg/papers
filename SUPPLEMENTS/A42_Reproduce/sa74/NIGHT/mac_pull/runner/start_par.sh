#!/bin/bash
# SA74: build and check the multi-core engine, then start runner 3 (plan A5: 20^4), which waits for runner 2's chains (about 15:20).
# Paste into a Terminal tab on the Mac, or ask your local Claude Code session to run it:
#     bash ~/repos/gr-qm-tests/sa74_night_regge/start_par.sh
cd "$(dirname "$0")"
clang -O2 -o regge_mac_par.new regge_r2_par.c -lpthread || { echo "build failed; nothing started"; exit 1; }
R2A=0.005 ./regge_mac_par.new run 3 0.2475 1 3 3 4242 chk_new.dat 0.5 - chk_new.cfg 10 0 0 > /dev/null
R2A=0.005 ./regge_mac_seq run 3 0.2475 1 3 3 4242 chk_old.dat 0.5 - chk_old.cfg 10 0 0 > /dev/null
if ! { cmp -s chk_new.cfg chk_old.cfg && diff <(grep -v '^#' chk_new.dat) <(grep -v '^#' chk_old.dat) > /dev/null; }; then
  echo "check 1 FAILED: without PAR the new engine differs from the current one; nothing started"; exit 1; fi
R2A=0.005 PAR=1 ./regge_mac_par.new run 12 0.248 1 1 1 777 chk_p1.dat 0.5 - chk_p1.cfg 10 0 0 > /dev/null
R2A=0.005 PAR=8 ./regge_mac_par.new run 12 0.248 1 1 1 777 chk_p8.dat 0.5 - chk_p8.cfg 10 0 0 > /dev/null
if ! cmp -s chk_p1.cfg chk_p8.cfg; then echo "check 2 FAILED: 1 and 8 threads give different runs; nothing started"; exit 1; fi
grep "edge order" chk_p8.dat
mv regge_mac_par.new regge_mac_par; rm -f chk_new.* chk_old.* chk_p1.* chk_p8.* STOP3
echo "engine checks passed: without PAR identical to the current engine; 1 and 8 threads identical"
pkill -f "[n]ight_queue3\.sh" 2>/dev/null; sleep 1
nohup caffeinate -is bash "$PWD/night_queue3.sh" >> night_queue3.log 2>&1 &
sleep 3; tail -n 2 night_queue3.log
echo "Runner 3 started: it waits until runner 2's chains finish (about 15:20), then runs plan A5 (20^4) with 8 threads. Gentle stop: touch ~/repos/gr-qm-tests/sa74_night_regge/STOP3"
