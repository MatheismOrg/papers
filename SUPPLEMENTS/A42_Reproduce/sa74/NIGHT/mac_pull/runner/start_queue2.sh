#!/bin/bash
# SA74 (Gary 03:24): keep the Mac computing past the night with the faster engine. Paste into a Terminal tab on the Mac:
#     bash ~/repos/gr-qm-tests/sa74_night_regge/start_queue2.sh
# Builds regge_mac_seq from regge_r2_big_x.c, checks that its default mode reproduces tonight's engine exactly, then starts runner 2
# in the background (kept awake). It takes over the cores as tonight's queue finishes. Gentle stop: touch the file STOP2 in this folder.
cd "$(dirname "$0")"
clang -O2 -o regge_mac_seq.new regge_r2_big_x.c || { echo "build failed; nothing started"; exit 1; }
R2A=0.005 ./regge_mac_seq.new run 3 0.2475 1 3 3 4242 chk_new.dat 0.5 - chk_new.cfg 10 0 0 > /dev/null
R2A=0.005 ./regge_mac run 3 0.2475 1 3 3 4242 chk_old.dat 0.5 - chk_old.cfg 10 0 0 > /dev/null
if cmp -s chk_new.cfg chk_old.cfg && diff <(grep -v '^#' chk_new.dat) <(grep -v '^#' chk_old.dat) > /dev/null; then
  mv regge_mac_seq.new regge_mac_seq; echo "engine check passed (default mode identical to tonight's engine)"
else echo "engine check FAILED; runner 2 not started"; exit 1; fi
rm -f chk_new.dat chk_new.cfg chk_old.dat chk_old.cfg STOP2
pkill -f "[n]ight_queue2\.sh" 2>/dev/null; sleep 1
touch queue2.txt
nohup caffeinate -is bash "$PWD/night_queue2.sh" >> night_queue2.log 2>&1 &
sleep 6; tail -n 3 night_queue2.log
echo "Runner 2 started: it takes over the cores as tonight's queue finishes, with no end time. Gentle stop: touch ~/repos/gr-qm-tests/sa74_night_regge/STOP2"
