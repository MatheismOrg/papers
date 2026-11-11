#!/bin/bash
# SA74 night: (re)start the server chains. Paste into a Terminal tab on the Mac:
#     bash ~/repos/gr-qm-tests/sa74_night_regge/start_server.sh
# Third version: the stop and the start now go in separate ssh calls. Before, the remote command line itself contained the
# script's path, so "pkill -f" matched its own shell and killed the start. The copy-back keeps running on the Mac.
cd "$(dirname "$0")"
KEY=~/.ssh/SERVER_KEY; HOST=root@SERVER_ADDRESS; SSH="ssh -i $KEY -o ConnectTimeout=15"
$SSH $HOST 'mkdir -p /root/sa74_night' || { echo "ssh failed"; exit 1; }
scp -q -i $KEY regge_x86 server_run.sh jobs_server.txt $HOST:/root/sa74_night/ || { echo "scp failed"; exit 1; }
$SSH $HOST 'pkill -f "[s]a74_night/server_run"; true'
$SSH $HOST 'cd /root/sa74_night && (nohup bash /root/sa74_night/server_run.sh > server_run.log 2>&1 < /dev/null &)'
sleep 5
$SSH $HOST 'cat /root/sa74_night/server_run.log; echo "engine processes running: $(pgrep -c regge_x86)"'
