#!/bin/bash
# SA74 night (Gary 00:11, server top): SA73's six DM-search workers also run at nice 19, so the six SA74 chains were taking
# about half of their CPU. This starts the Mac queue runner (the big grids move to the Mac's cores as they free up), stops the
# SA74 chains on the server, and stops the Mac copy-back. Paste into a Terminal tab on the Mac (not inside ssh):
#     bash ~/repos/gr-qm-tests/sa74_night_regge/switch_to_mac.sh
cd "$(dirname "$0")"
KEY=~/.ssh/SERVER_KEY; HOST=root@SERVER_ADDRESS; SSH="ssh -i $KEY -o ConnectTimeout=15"
bash start_queue.sh
$SSH $HOST 'P=$(pgrep -f "[s]a74_night/server_run" | tr "\n" "," | sed "s/,$//"); [ -n "$P" ] && { pkill -P "$P"; kill $(echo $P | tr "," " "); }; for i in 1 2 3; do pkill -f "[s]a74_night/regge_x86"; sleep 2; done; echo "server: SA74 engine processes left: $(pgrep -f "[s]a74_night/regge_x86" | wc -l)"; nvidia-smi -L 2>/dev/null || echo "server: no NVIDIA GPU found"' || echo "ssh to the server failed"
pkill -f "[s]erver/copyback.sh" && echo "Mac copy-back stopped."
