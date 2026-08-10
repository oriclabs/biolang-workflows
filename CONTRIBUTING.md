# Contributing

Keep workflows practical, reproducible, and explicit about their scientific
contract.

- Pin BioLang/package versions and random seeds.
- Record data source, accession, licence, retrieval date, and SHA-256 hashes.
- Do not commit large downloaded datasets or generated validation results.
- Keep external R/Python tools optional and out of the BioLang process.
- Label oracle observations separately from BioLang implementation code.
- Do not translate or copy copyleft implementation source into MIT BioLang code.
- State CPU/GPU backend and thread count when results can differ.
- Prefer a small smoke fixture for pull requests and a scheduled real-data run.
- Do not claim scientific equivalence from correlation alone; include
  scale-sensitive error, ranking, and downstream metrics where appropriate.

Run `python scripts/check_repository.py` before submitting a change.

