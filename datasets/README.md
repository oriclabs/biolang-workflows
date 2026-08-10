# Dataset policy

This directory contains metadata and download instructions, not large datasets.
Each dataset manifest should include:

- stable source URL or accession;
- upstream authors and citation;
- data licence and redistribution status;
- retrieval date;
- SHA-256 for every downloaded input;
- expected matrix dimensions and format;
- deterministic filtering or preparation commands;
- whether the data contain human or otherwise sensitive information.

Downloaded files belong under `datasets/cache/` or a user-selected external
directory; both must remain untracked. A workflow must fail clearly when a hash
or expected dimension does not match.

