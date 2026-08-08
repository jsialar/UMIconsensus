# TODO

## Handle samples that genuinely carry an indel at a SNP position

The SNP sites in `--snplist` were selected on the assumption that none of them
carry an indel, based on checking the positions against the HGDP+1KG genomic
data. `call_consensus.py` therefore treats a `D` or `I` at any of those
positions as evidence of an ONT error or misalignment and discards the whole
UMI family.

A rare or previously unobserved sample could genuinely carry an indel at one of
these sites. Today the pipeline would silently drop the families carrying it —
the true variant never reaches `_haps.csv`, and nothing in the output signals
that it happened.

Possible directions:

- Report, rather than only filter: emit a per-target count of families dropped
  for `D`/`I`, broken down by position, so a consistent indel across many
  independent families is visible without re-running.
- Decide what a real indel should look like in the output. The minihaplotype
  strings are currently one character per SNP and directly comparable; letting
  an indel through would either need a code that preserves that alignment or a
  representation that tolerates variable-length haplotypes.
  (`process_consensus.py` already sizes its quality rows per haplotype rather
  than assuming a shared length, so it is partly ready for this.)
- Consider whether a family with an indel at one site should still contribute
  its calls at the other sites, instead of being discarded entirely.
