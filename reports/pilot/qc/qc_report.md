# QC report

Genomes requested: **2000** (Bacteria: 1796, Archaea: 204)

| step | Archaea | Bacteria | all | % of requested |
|---|---|---|---|---|
| fetched | 204 | 1796 | 2000 | 100.0 |
| has_16s | 204 | 1795 | 1999 | 100.0 |
| annotated | 204 | 1796 | 2000 | 100.0 |
| in_features | 204 | 1795 | 1999 | 100.0 |
| usable | 204 | 1795 | 1999 | 100.0 |
| pcr_V3V4 | 158 | 1759 | 1917 | 95.9 |
| pcr_V4 | 185 | 1776 | 1961 | 98.0 |

`usable` = in the 16S feature matrix and has at least one KO annotation; these genomes enter modelling.
`pcr_<arm>` = in-silico PCR success with that primer pair.

Design expectation: usable >= 90 % of requested. Observed: **100.0 %**.

## Not usable (1)

| reason | n |
|---|---|
| no full-length 16S | 1 |
