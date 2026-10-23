# A42 supplement: the rule search

These files reproduce the search in "Searching the Shortest Rewriting Rules for Four-Dimensional Spacetime: A First Pass and an Obstacle". Python 3.10 or later, standard library only.

## Files

- A42_Rule_Search.py: the rewriting engine, the causal graph, the three dimension tests and the coarse-grained check, the validation on sprinkled flat intervals, the link counts and the two grid controls.
- A42_Compile_Results.py: compiles the run files into the two data files below.
- A42_Search_Results.json: the validation, the census of all four classes of rules, the 27 long runs of the first pass and the ten of the second with their size series, and the de Sitter comparison.
- A42_All_Rules.csv: one row per rule, 6,435 rows: its class, length in bits, behavior, updates, generations, growth fits, and the first two tests on its largest and inner intervals.

## How to reproduce

    python3 A42_Rule_Search.py validate validation.json
    python3 A42_Rule_Search.py enumerate rules_22_32.json
    python3 A42_Rule_Search.py run rules_22_32.json res_all.jsonl --events 3000
    python3 A42_Rule_Search.py sample rules_22_42_sample.json 1500
    python3 A42_Rule_Search.py run rules_22_42_sample.json s42_all.jsonl --events 3000
    python3 A42_Rule_Search.py run deep2_rules.json deep2.jsonl --events 20000 --sec 60 --series
    python3 A42_Rule_Search.py run timelimited_rules.json tl_a.jsonl --events 3000 --sec 30
    python3 A42_Rule_Search.py sample2 rules_32_42_sample.json 3 4 2 1000 42
    python3 A42_Rule_Search.py run2 rules_32_42_sample.json s3242.jsonl --events 3000
    python3 A42_Rule_Search.py sample2 rules_23_33_sample.json 2 3 3 1000 42
    python3 A42_Rule_Search.py run2 rules_23_33_sample.json s2333.jsonl --events 3000
    python3 A42_Rule_Search.py run2 timelimited2_rules.json tl2_a.jsonl --events 3000 --sec 30
    python3 A42_Rule_Search.py long long_rules.json long.jsonl --events 120000 --sec 600
    python3 A42_Rule_Search.py desitter desitter.json
    python3 A42_Compile_Results.py OUT_DIR

deep2_rules.json holds the 27 rules listed under long_runs in A42_Search_Results.json, each with its group, and long_rules.json the ten listed under long_runs_second_pass. The second pass's commands (sample2, run2, long) use the general engine, run_general, which handles relations of any size and left sides of any length and reproduces run exactly on two-relation left sides; a long run reads its tests inside the past of a top update whose past holds at most 70,000 updates, so the stored pasts fit in memory. timelimited2_rules.json holds the 18 second-pass rules that hit the 8-second limit below 300 updates. timelimited_rules.json holds the 12 slow single chains that hit the 8-second limit below 300 updates; the compile step counts them by their 30-second reruns. The run command takes --start and --stop to split a list across processors.

## Notes

Results are deterministic apart from the time limits. The scan rows were made by an earlier version of the script, before the longest-chain test and the coarse-grained check were added; the engine and the interval choices are unchanged, and on a random 60 rules of the full class the current script reproduces 59 rows exactly (the other is one of the 12 slow chains). The field cone_dim, the growth of the fullest past with depth, reads low on sprinkled four-dimensional intervals at these sizes (about 2.4 to 2.6), so the paper does not use it.
