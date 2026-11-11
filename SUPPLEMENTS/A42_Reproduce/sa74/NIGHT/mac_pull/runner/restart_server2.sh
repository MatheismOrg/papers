#!/bin/bash
# SA74 plan S2, restart (12:00 PST). The first start (11:57) failed at once: all six engines got "Text file busy", because the freshly
# copied regge_x86 was still held open on the server. Nothing ran, so nothing is lost. This gives the engine a fresh file and starts the
# same six jobs again, at idle priority, if the server is still mostly idle. Paste into a Terminal tab on the Mac:
#     bash ~/repos/gr-qm-tests/sa74_night_regge/restart_server2.sh
cd "$(dirname "$0")"
KEY=~/.ssh/SERVER_KEY; HOST=root@SERVER_ADDRESS; SSH="ssh -i $KEY -o ConnectTimeout=15"
read LOAD NCPU < <($SSH $HOST 'echo $(cut -d" " -f1 /proc/loadavg) $(nproc)') || { echo "ssh failed"; exit 1; }
echo "server: load $LOAD on $NCPU cores"
if ! awk -v l="$LOAD" -v n="$NCPU" 'BEGIN { exit !(l < n / 2) }'; then echo "server busy (SA73's search?): not started"; exit 0; fi
$SSH $HOST 'cd /root/sa74_s2 && cp regge_x86 regge_x86.new && mv -f regge_x86.new regge_x86 && chmod +x regge_x86' || { echo "engine refresh failed"; exit 1; }
$SSH $HOST 'cd /root/sa74_s2 && (P=6 nohup chrt -i 0 bash /root/sa74_s2/server_run.sh > server_run.log 2>&1 < /dev/null &)'
sleep 10
$SSH $HOST 'cat /root/sa74_s2/server_run.log; echo "engines running: $(pgrep -fc "[s]a74_s2/regge_x86")"'
