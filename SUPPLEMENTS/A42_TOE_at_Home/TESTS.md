# TOE@home: tests before release (November 11, 2026)

All tests ran on Linux with Python 3.11, 3.12 and 3.13.

## The computation

1. **Same rules as the paper.** The blocks of class b2_3 list exactly the 2,935 rules of the paper's complete class, each
   once, in the form the paper runs.
2. **Same results as the paper.** A full run of class b2_3 (2,935 rules, 657 seconds on one core) equals the paper's
   stored record for 2,905 rules on all eight fields compared (updates, generations, and size, ordering fraction
   dimension and midpoint dimension of both intervals). The other 30 are the rules cut off by the step budget here and
   by the 8-second clock in the paper; all 30 grow chains or trees.
3. **Long runs as in the paper.** Four of the paper's long runs (20,000 updates) were repeated; every field and the whole
   size series equal the paper's stored record.
4. **Same result on every Python version.** One block gave the same digest on Python 3.11, 3.12 and 3.13.
5. **Indexing.** The direct indexing of right sides gives the order of itertools.combinations_with_replacement, and the
   direct list of left sides equals the paper program's list for three classes.

## The coordinator, without the network (30 checks, all passed)

Real block results, several pretend volunteers, one of them faulty:

- two volunteers at one address never check each other; a volunteer elsewhere gets the same block
- a wrong token, a missing rule, another rule's text, a wrong field type and another version are refused
- unknown fields are dropped before a result is stored
- a different reading, and a count off by one, become a dispute; a third volunteer settles it; the faulty volunteer is
  marked and gets no credit
- readings 0.001 apart agree
- an expired reservation is handed out again
- a trusted machine validates alone
- a near reading creates a long run of that rule, handed out ahead of everything else and checked the same way (tested
  with the window moved to two dimensions, where the paper's class has rules)
- an offline result file is taken in; the progress page renders; exported records equal the trusted run

## Over the network on one machine (18 checks, all passed)

- two volunteers on Python 3.13 and 3.11 check five blocks; a third is told the work is finished
- a nickname with markup is cleaned before it reaches the page
- garbage, an unknown address, an oversize upload, a compressed bomb, another version and a made-up result are refused,
  and the coordinator stays up
- the self-test passes
- a client killed in mid-block keeps its finished rules and completes the block after a restart
- two processes on one machine work side by side
- an offline run keeps its result file, uploads it later, and a second upload of the same file adds nothing

## Not tested

- Windows and macOS (the program uses nothing specific to an operating system, apart from lowering its priority)
- many volunteers at once; a public address; HTTPS, which the web server in front provides
