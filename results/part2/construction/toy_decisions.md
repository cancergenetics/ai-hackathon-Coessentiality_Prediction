| Variant | Planted complexes recovered (F1 >= 0.8, any level; planted + hier) | of which size 3-4 | Sub-modules recovered | hier parents at level 1 | hier overlapping complexes | False clusters (< 50% from one complex) | Null-scenario clusters | Arm-block non-locus clusters |
|---|---|---|---|---|---|---|---|---|
| V0_jaccard_greedy | 67/105 | 6/30 | 39/39 | 1/9 | 12/12 | 0 | 0 | 0 of 19 |
| V1_containment_greedy | 85/105 | 24/30 | 39/39 | 2/9 | 12/12 | 6 | 0 | 1 of 25 |
| V2_containment_nested | 86/105 | 25/30 | 39/39 | 2/9 | 12/12 | 6 | 0 | 1 of 25 |
| V3_jaccard_nested | 67/105 | 6/30 | 39/39 | 1/9 | 12/12 | 0 | 0 | 0 of 19 |

Attachment (D3, V2 map):

| Toy seed | Scenario | Held-out members | Attached to own complex | Wrong-cluster pairs | Other genes attached (by class) | Tests |
|---|---|---|---|---|---|---|
| 0 | planted | 10 | 10 | 1 | {'never': 1} | 30702 |
| 0 | hier | 17 | 17 | 0 | {'never': 2} | 61642 |
| 0 | null, 14 random size-matched sets | 0 | - | - | 0 attached | 25144 |
| 0 | armblock | 0 | 0 | 0 | {'never': 92} | 12572 |
| 1 | planted | 10 | 10 | 2 | {'never': 1} | 30668 |
| 1 | hier | 17 | 17 | 0 | {'never': 2} | 57952 |
| 1 | null, 14 random size-matched sets | 0 | - | - | 0 attached | 25116 |
| 1 | armblock | 0 | 0 | 0 | {'never': 114} | 14352 |
| 2 | planted | 10 | 10 | 0 | {} | 28896 |
| 2 | hier | 17 | 17 | 0 | {} | 58016 |
| 2 | null, 14 random size-matched sets | 0 | - | - | 0 attached | 25144 |
| 2 | armblock | 0 | 0 | 0 | {'never': 116} | 17960 |
