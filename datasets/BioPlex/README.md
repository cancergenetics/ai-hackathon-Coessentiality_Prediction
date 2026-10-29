# BioPlex 3.0 (placeholder)

BioPlex 3.0 interaction networks (Huttlin et al. 2021, *Cell*), <https://bioplex.hms.harvard.edu/>. Used only for the pre-registered BioPlex bridge (validation).

Put these files here. They are not distributed with this repository and are git-ignored; do not commit them.

| File | Version | SHA-256 of the file used | Licence (as recorded in DATA.md) |
|---|---|---|---|
| `BioPlex_293T_Network_10K_Dec_2019.tsv` | Dec 2019 release | `5a98def087fc77c7bf07aac79e2553fd64ea1a27a7d0c4617870476f3f120833` | BioPlex terms (free for research; unverified) |
| `BioPlex_HCT116_Network_5.5K_Dec_2019.tsv` | Dec 2019 release | `b37be6de6fe15e948151ab9480cfcb88a6795780cf7a4059b386b07f3c1e760d` | as above |

Check with `shasum -a 256 <file>` (or `scripts/check_data.sh`). Other releases will give somewhat different results.
