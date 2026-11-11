#!/bin/bash
# SA74 night (Gary 23:11): run the big-grid gravity chains natively on this Mac, on every core but one.
# Paste into Terminal:   bash ~/repos/gr-qm-tests/sa74_night_regge/start_native.sh
# It starts in the background, kept awake with caffeinate, and returns; you can close the window. Running it again stops the
# running copy and resumes from the saved configurations. Progress: ~/repos/gr-qm-tests/sa74_night_regge/night_native.log
# Stop everything:   pkill -f sa74_night_regge
cd "$(dirname "$0")"
command -v clang >/dev/null || { echo "No clang here. Install it with: xcode-select --install"; exit 1; }
pkill -f "$PWD/night_native.sh" 2>/dev/null; pkill -f "$PWD/regge_mac" 2>/dev/null; sleep 1
nohup caffeinate -is bash "$PWD/night_native.sh" > night_native.log 2>&1 &
echo "Started on $(sysctl -n hw.ncpu) cores (process $!). Progress: ~/repos/gr-qm-tests/sa74_night_regge/night_native.log. You can close this window."
