#!/bin/bash
# SA74 plan S2, fixed 11:50 PST before its runs (Gary 11:40: "how can it finish faster? run on both mac and server?"). Paste into a
# Terminal tab on the Mac (not inside ssh):
#     bash ~/repos/gr-qm-tests/sa74_night_regge/start_server2.sh
# Six more forks of plan A2's two 16^4 chains at Hamber's k_c = 0.244 (three from each), on the server, each 500 sweeps at the coupling
# and then 1,500 measured, with the static engine regge_x86. It starts only if the server is mostly idle (load below half its cores),
# and it runs at idle priority (chrt -i 0), so SA73's search keeps its CPU. A copy-back brings the results into server2/ here every
# 10 minutes until 22:00. Stop on the server: ssh -i ~/.ssh/SERVER_KEY root@SERVER_ADDRESS 'pkill -f "[s]a74_s2"'
cd "$(dirname "$0")"
KEY=~/.ssh/SERVER_KEY; HOST=root@SERVER_ADDRESS; SSH="ssh -i $KEY -o ConnectTimeout=15"
read LOAD NCPU < <($SSH $HOST 'echo $(cut -d" " -f1 /proc/loadavg) $(nproc)') || { echo "ssh failed"; exit 1; }
echo "server: load $LOAD on $NCPU cores"
if ! awk -v l="$LOAD" -v n="$NCPU" 'BEGIN { exit !(l < n / 2) }'; then echo "server busy (SA73's search?): not started"; exit 0; fi
mkdir -p server2
for i in 1 2 3; do for p in r1:a r2:b; do                      # whole copies of the parents' latest configurations
  src=A2_L16_k0.244_${p%%:*}.cfg; dst=server2/S2_L16_k0.244_${p##*:}$i.cfg
  for t in 1 2 3; do cp $src $dst; cmp -s $src $dst && [ $(stat -f %z $dst) = 7864332 ] && break; sleep 20; done
done; done
cat > server2/jobs_server.txt <<'EOJ'
# SA74 plan S2 (see start_server2.sh): TAG L K THERM MEAS SEED A
S2_L16_k0.244_a1 16 0.244 500 1500 99930 0.005
S2_L16_k0.244_b1 16 0.244 500 1500 99931 0.005
S2_L16_k0.244_a2 16 0.244 500 1500 99932 0.005
S2_L16_k0.244_b2 16 0.244 500 1500 99933 0.005
S2_L16_k0.244_a3 16 0.244 500 1500 99934 0.005
S2_L16_k0.244_b3 16 0.244 500 1500 99935 0.005
EOJ
$SSH $HOST 'mkdir -p /root/sa74_s2' && scp -q -i $KEY regge_x86 server_run.sh server2/jobs_server.txt server2/S2_*.cfg $HOST:/root/sa74_s2/ || { echo "scp failed"; exit 1; }
$SSH $HOST 'cd /root/sa74_s2 && (P=6 nohup chrt -i 0 bash /root/sa74_s2/server_run.sh > server_run.log 2>&1 < /dev/null &)'
sleep 8
$SSH $HOST 'cat /root/sa74_s2/server_run.log; echo "engines running: $(pgrep -fc "[s]a74_s2/regge_x86")"'
cat > server2/copyback.sh <<EOS
#!/bin/bash
# plan S2 copy-back: every 10 minutes until 22:00, the server's chain files (not the per-segment parts) and its log, into server2/
while [ \$(date +%H%M) -lt 2200 ]; do
  scp -q -i $KEY "$HOST:/root/sa74_s2/S2_L16_k0.244_[ab][123].dat" "$HOST:/root/sa74_s2/S2_L16_k0.244_[ab][123].log" "$HOST:/root/sa74_s2/server_run.log" "$PWD/server2/" 2>/dev/null
  sleep 600
done
EOS
nohup caffeinate -is bash "$PWD/server2/copyback.sh" > server2/copyback.log 2>&1 &
echo "Copy-back running: the server's results land in ~/repos/gr-qm-tests/sa74_night_regge/server2 every 10 minutes."
