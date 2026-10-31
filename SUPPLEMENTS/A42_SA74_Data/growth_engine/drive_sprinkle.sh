cd /home/claude/sa74/big
for s in 1 2 3; do ./tp_big sprinkle 4 8000 $s; done >> sprinkle4.jsonl
for s in 1 2; do ./tp_big sprinkle 4 32000 $s; done >> sprinkle4.jsonl
for s in 1 2; do ./tp_big sprinkle 4 64000 $s; done >> sprinkle4.jsonl
./tp_big sprinkle 4 128000 1 >> sprinkle4.jsonl
echo DONE >> sprinkle4.jsonl
