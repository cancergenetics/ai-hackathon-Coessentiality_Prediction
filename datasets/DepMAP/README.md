# DepMap (placeholder)

DepMap Public 26Q3 files, <https://depmap.org/portal/data_page/>. `CRISPRGeneEffect.csv` is the outcome and the input of the locked method; the Omics files, `Model.csv`, `ModelCondition.csv` and `CRISPRConfounders.csv` define contexts and covariates; `ScreenGeneEffect.csv` is used only for the screen-quality sensitivity.

Put these files here. They are not distributed with this repository and are git-ignored; do not commit them.

| File | Version | SHA-256 of the file used | Licence (as recorded in DATA.md) |
|---|---|---|---|
| `CRISPRGeneEffect.csv` | 26Q3 (as stated by user; files carry no version stamp) | `1e0c492944a7297c0eda179081a1e81a1cfe41a97628facd7eddbdda1a7504ee` | CC BY 4.0 (DepMap public data terms) |
| `ScreenGeneEffect.csv` | 26Q3 (user) | `08e15ab81d7450db70c48c7781f312f6cebc068da7077bea775569fb516a3f04` | CC BY 4.0 |
| `CRISPRConfounders.csv` | 26Q3 (user) | `0e098325d53403d7852f3a773477d20a508a9bf29f43e10fe28bf0f550b35779` | CC BY 4.0 |
| `Model.csv` | 26Q3 (user) | `65eb1895800e71a963a3ab56e5d08c8d006cb0d9e7aa68f7e34a3810a2a3cbc7` | CC BY 4.0 |
| `ModelCondition.csv` | 26Q3 (user) | `a1edda78cd60718e5c99690ce8ab80315b65ddc9346643338d2c1c774fc51326` | CC BY 4.0 |
| `OmicsSomaticMutationsMatrixDamaging.csv` | 26Q3 (user-stated) | `736b7a82b4d1f234dc597a409d104186c282f01b0caa56365c27ff987ae49e20` | CC BY 4.0 |
| `OmicsSomaticMutationsMatrixHotspot.csv` | 26Q3 (user) | `0cbdacdce2355d2dd600d3741a36ae093dbe77d7bee44bde97411b1e96973058` | CC BY 4.0 |
| `OmicsCNGeneWGS.csv` | 26Q3 (user); ploidy-relative linear CN ratio (diploid ≈ 1), not log2 (round-1 reviews) | `31cace66f202adaace306d975b5661c8b77ea82f07d2ed34045eec9c7958fcd6` | CC BY 4.0 |
| `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | 26Q3 (user) | `252d1fedc7014953b19ffe5d68151b4df973c8078b12c9e5ccd97e292548d177` | CC BY 4.0 |
| `OmicsProfiles.csv` | 26Q3 (user) | `4241e537f0b23139c0e525b4d78fb53ee65c5507771a33e2cd5debd78f7179ea` | CC BY 4.0 |
| `OmicsGlobalSignatures.csv` | 26Q3 (user) | `51d242e99176e82238df3c90a62b848b40a459a4a4981a2b0a923ed46ad8fb93` | CC BY 4.0 |
| `OmicsFusionFiltered.csv` | 26Q3 (user) | `985ac6383a95c39f36856c583f32c7990c9f5464d982a2dd84a9dc07fcbc979b` | CC BY 4.0 |

Check with `shasum -a 256 <file>` (or `scripts/check_data.sh`). Other releases will give somewhat different results.
