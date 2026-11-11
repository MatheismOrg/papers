#!/bin/bash
# SA74 plan S3 (fixed 13:05 PST): a second 20^4 chain at k = 0.244 on the server with the multi-core engine, 6 threads at idle priority
# (the other project and SA73 keep their CPU; nothing in /root/OTHER_PROJECT is touched). It starts only if the server's load is below half its cores
# and at least 4 GB of memory are free. A copy-back brings the results into server3/ here every 10 minutes until midnight.
# Paste into a Terminal tab on the Mac (or ask your local Claude Code session):   bash ~/repos/gr-qm-tests/sa74_night_regge/start_server3.sh
# Stop on the server: ssh -i ~/.ssh/SERVER_KEY root@SERVER_ADDRESS 'pkill -f "[s]a74_s3"'
cd "$(dirname "$0")"
KEY=~/.ssh/SERVER_KEY; HOST=root@SERVER_ADDRESS; SSH="ssh -i $KEY -o ConnectTimeout=15"
read LOAD NCPU FREE < <($SSH $HOST 'echo $(cut -d" " -f1 /proc/loadavg) $(nproc) $(free -m | awk "/^Mem:/ {print \$7}")') || { echo "ssh failed"; exit 1; }
echo "server: load $LOAD on $NCPU cores, $FREE MB free"
if ! awk -v l="$LOAD" -v n="$NCPU" -v f="$FREE" 'BEGIN { exit !(l < n / 2 && f > 4000) }'; then echo "server busy or short of memory: not started"; exit 0; fi
$SSH $HOST 'mkdir -p /root/sa74_s3' && scp -q -i $KEY regge_par_x86 $HOST:/root/sa74_s3/regge_par_x86.up && scp -q -i $KEY server_run_par.sh jobs_s3.txt $HOST:/root/sa74_s3/ || { echo "scp failed"; exit 1; }
$SSH $HOST 'cd /root/sa74_s3 && cp regge_par_x86.up regge_par_x86 && chmod +x regge_par_x86 && (nohup chrt -i 0 bash /root/sa74_s3/server_run_par.sh > server_run.log 2>&1 < /dev/null &)'
sleep 15
$SSH $HOST 'cat /root/sa74_s3/server_run.log; echo "engine threads running: $(ps -eLf | grep -c "[s]a74_s3/regge_par_x86")"; grep -q "^SA74 .*sa74_s3" /root/RUNNING.md 2>/dev/null || echo "SA74 | MATHEISM WORKS (A42 lattice gravity) | plan S3: one 20^4 Regge chain, /root/sa74_s3 (multi-core engine) | 6 threads | idle (chrt -i 0) | $(date +"%b %-d %H:%M") server time | about 5 hours | pkill -f \"[s]a74_s3\"" >> /root/RUNNING.md'
mkdir -p server3
cat > server3/copyback.sh <<EOS
#!/bin/bash
# plan S3 copy-back: every 10 minutes until midnight, the chain file, its log and the runner log, into server3/
while [ \$(date +%H%M) -lt 2359 ]; do
  scp -q -i $KEY "$HOST:/root/sa74_s3/S3_L20_k0.244_r2.dat" "$HOST:/root/sa74_s3/S3_L20_k0.244_r2.log" "$HOST:/root/sa74_s3/server_run.log" "$PWD/server3/" 2>/dev/null
  sleep 600
done
EOS
nohup caffeinate -is bash "$PWD/server3/copyback.sh" > server3/copyback.log 2>&1 &
echo "Copy-back running: results land in ~/repos/gr-qm-tests/sa74_night_regge/server3 every 10 minutes."
