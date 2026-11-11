#!/bin/bash
# Run logs for the shared machines (Gary 12:05, Nov 7): ~/RUNNING/RUNNING.md on this Mac and /root/RUNNING.md on the server,
# each started with a header and SA74's lines (other sessions add theirs). If a log already exists, only SA74's lines are appended.
# Paste into a Terminal tab on the Mac:   bash ~/repos/gr-qm-tests/sa74_night_regge/make_run_logs.sh
cd "$(dirname "$0")"
KEY=~/.ssh/SERVER_KEY; HOST=root@SERVER_ADDRESS; SSH="ssh -i $KEY -o ConnectTimeout=15"
mkdir -p ~/RUNNING
if [ -f ~/RUNNING/RUNNING.md ]; then tail -n +5 run_log_mac.md >> ~/RUNNING/RUNNING.md; else cp run_log_mac.md ~/RUNNING/RUNNING.md; fi
if $SSH $HOST 'test -f /root/RUNNING.md'; then
  scp -q -i $KEY run_log_server.md $HOST:/root/RUNNING.sa74.tmp && $SSH $HOST 'tail -n +5 /root/RUNNING.sa74.tmp >> /root/RUNNING.md; rm -f /root/RUNNING.sa74.tmp'
else
  scp -q -i $KEY run_log_server.md $HOST:/root/RUNNING.md
fi
echo "== Mac: ~/RUNNING/RUNNING.md"; cat ~/RUNNING/RUNNING.md
echo; echo "== Server: /root/RUNNING.md"; $SSH $HOST 'cat /root/RUNNING.md; echo; uptime'
