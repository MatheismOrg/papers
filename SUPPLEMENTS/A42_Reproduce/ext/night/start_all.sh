#!/bin/bash
# SA74 night (Gary 23:24): one paste starts everything. Paste into a NEW Terminal tab on the Mac (Cmd+T), not inside ssh:
#     bash ~/repos/gr-qm-tests/sa74_night_regge/start_all.sh
# 1. starts the Mac chains natively (start_native.sh); 2. copies the static engine and plan S1 to the server and starts them
# there at low priority; 3. starts a background copy-back that brings the server's results into server/ here every 10 minutes
# until 06:00. You can close the tab afterwards. Stop everything: pkill -f sa74_night (here) and the same on the server.
cd "$(dirname "$0")"
KEY=~/.ssh/SERVER_KEY; HOST=root@SERVER_ADDRESS; SSH="ssh -i $KEY -o BatchMode=yes -o ConnectTimeout=15"
bash start_native.sh
$SSH $HOST 'mkdir -p ~/sa74_night' && scp -q -i $KEY regge_x86 server_run.sh jobs_server.txt $HOST:sa74_night/ \
  && $SSH $HOST 'cd ~/sa74_night && pkill -f "sa74_night/server_run.sh"; nohup bash ~/sa74_night/server_run.sh > server_run.log 2>&1 < /dev/null &' \
  && echo "Server started. Its log: ~/sa74_night/server_run.log on the server." || echo "Server step failed (ssh or scp)."
mkdir -p server
cat > server/copyback.sh <<EOS
#!/bin/bash
# every 10 minutes until 06:00: the server's chain files (not the per-segment parts) and logs, into this folder
cd "\$(dirname "\$0")"
while [ "\$(date +%H)" != 06 ]; do
  $SSH $HOST 'cd ~/sa74_night && tar czf - \$(ls S_*.dat | grep -v part) *.log 2>/dev/null' | tar xzf - 2>/dev/null && echo "\$(date +%H:%M) copied"
  sleep 600
done
EOS
nohup caffeinate -is bash "$PWD/server/copyback.sh" > server/copyback.log 2>&1 &
echo "Copy-back running (process $!): server results land in ~/repos/gr-qm-tests/sa74_night_regge/server every 10 minutes."
