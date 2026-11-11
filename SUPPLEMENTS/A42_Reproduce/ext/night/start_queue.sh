#!/bin/bash
# SA74 night (Gary 00:06): keep every Mac core busy all night. Paste into a Terminal tab on the Mac:
#     bash ~/repos/gr-qm-tests/sa74_night_regge/start_queue.sh
# Starts the queue runner in the background (kept awake with caffeinate) and returns; the window can be closed.
cd "$(dirname "$0")"
pkill -f "[n]ight_queue.sh" 2>/dev/null; sleep 1
touch queue.txt
nohup caffeinate -is bash "$PWD/night_queue.sh" >> night_queue.log 2>&1 &
sleep 12; tail -n 4 night_queue.log
echo "Mac queue runner started: it starts the queued chains as cores free up."
