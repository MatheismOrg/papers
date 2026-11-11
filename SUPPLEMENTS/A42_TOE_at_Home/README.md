# TOE@home

A volunteer search for a short rule that grows four-dimensional spacetime. It continues the search of the paper
"Shortest Rules Search" (A42) on the computers of anyone who wants to help, the way SETI@home shared its work.

## What it looks for

A rule rewrites a few relations between points into a few new ones, again and again. Every rewrite is an update. The
order of cause and effect among the updates is a candidate spacetime. The program runs each rule for up to 3,000
updates and reads the dimension of what grew with three tests from causal-set theory: the ordering fraction, the
midpoint and the longest chain. Our spacetime reads four on all three.

The paper tested all 2,935 rules of the smallest class and samples of three larger ones. None reads four on all three
tests. This project tests the larger classes completely.

## The odds, plainly

Low. The paper gives a reason why rules of this kind struggle: each update has a limited number of direct causes, and
four dimensions need ever more. A hit would be a major result. The likely outcome is a wider negative result: no rule
in these classes grows four dimensions. That is worth having too, because it is complete and checked twice.

## For volunteers

You need Python 3.10 or later. Nothing else is installed.

1. Unpack the folder. Keep `toehome.py` and `A42_Rule_Search.py` together.
2. Run `python3 toehome.py` (on Windows: `py toehome.py`).
3. Stop with Ctrl-C at any time. The next start continues where it stopped.

Options:

- `--cores 4` uses four processes.
- `--nick Ada` puts a label on the list of contributors.
- `--selftest` checks in about a minute that your machine reproduces a known block.
- `--offline` works without a network and keeps result files in `toehome_data/done`; `--send` uploads them later.

What the program does on your machine:

- It runs at the lowest priority, so other programs come first.
- It contacts one address, the coordinator, to fetch the label of a block (for example `b2_4:3:1:3000`) and to send
  the result. It sends a random label for your machine, your nickname if you gave one, and the result. Nothing else.
- It cannot receive code. A block label only selects which rules to run; the computation is fixed in the two files.
- It checks the SHA-256 of `A42_Rule_Search.py` before use. That file is the search program published with the paper,
  unchanged.
- It writes only inside the folder `toehome_data` beside the program.

## How results are checked

- **Same result on every machine.** There is no time limit inside a run, only a fixed budget of matching steps. Tested
  on Python 3.11, 3.12 and 3.13.
- **Every block is done twice,** by volunteers at different network addresses. Whole numbers must be equal and
  readings equal to 0.002. A third volunteer settles a disagreement.
- **The coordinator works out the rules of each block itself** and refuses a result that lists anything else.
- **Two stages.** A rule that reads between 3.5 and 4.5 on all three tests gets a long run of 20,000 updates, checked
  twice the same way. The long run is judged by the criterion the paper fixed in advance: the largest interval of at
  least 1,000 elements and the two largest of the size series must each read four within three standard deviations.
- **A pass is not a claim.** Any rule that passes must be rerun by hand, at larger sizes and under other update
  orders, before anything is said about it.

## The classes

| Class | Left side | Right side | Rules | One core |
|---|---|---|---|---|
| b2_3 | 2 relations | 3 relations | 2,935 | 17 minutes |
| b2_4 | 2 relations | 4 relations | 18,263 | 2.4 hours |
| b2_5 | 2 relations | 5 relations | about 108,000 | 7 hours |
| b3_4 | 3 relations | 4 relations | about 147,000 | 12 hours |
| b3_5 | 3 relations | 5 relations | about 1.8 million | 23 days |
| t2_3 | 2 three-node relations | 3 three-node relations | about 42 million | 100 days |
| b4_5 | 4 relations | 5 relations | about 30 million | 230 days |
| t2_4 | 2 three-node relations | 4 three-node relations | about 3 billion | 60 years |

Rule counts from the third row on are estimates from samples. Times are for one core, measured on a cloud server, and
double when every block is done twice. Class b2_3 is the paper's own class: it is kept as the first job so that the project
starts by reproducing a published result.

## For whoever runs the coordinator

The coordinator is `toehome_server.py`. It needs Python 3.10 and nothing else. Its data is one SQLite file and one
folder of result files.

    python3 toehome_server.py init  --db toe.sqlite
    python3 toehome_server.py serve --db toe.sqlite --port 8742
    python3 toehome_server.py report --db toe.sqlite

- By default it listens on 127.0.0.1 only. Put it behind a web server with HTTPS and start it with `--behind-proxy`.
- Put the public address into `toehome.json` beside `toehome.py` before you publish the package:
  `{"server": "https://your.address"}`.
- `trust CLIENT_ID` marks a machine of your own; its results count without a second one. The label of a machine is in
  its `toehome_data/client.json`.
- `load FILE...` takes in result files from volunteers who ran offline.
- `export OUT.jsonl` writes every checked rule record.
- It stores no network address. It keeps a salted hash of the address prefix, used only to make sure the two results
  of a block come from different places.

## What is not in this version

The gravity lattice and the CDT runs of the paper are compiled programs with random sampling. Sharing them needs builds
for Windows, Mac and Linux and a statistical check of results, since two machines do not produce the same chain. That
is a separate piece of work.

## Files

- `toehome.py`: the volunteer program.
- `A42_Rule_Search.py`: the paper's search program, unchanged
  (SHA-256 4739063aceab87cecb1054c12c603c2ea9d433a2773962bfcbc7ab05a4a2418e).
- `toehome_server.py`: the coordinator.
- `TESTS.md`: what was tested before release, with results.
