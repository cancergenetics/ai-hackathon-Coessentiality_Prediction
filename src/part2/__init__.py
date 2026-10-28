"""Part 2 complex-construction infrastructure (DESIGN.md v2 §2.1, §2.3).

Import as the top-level package ``part2`` (``src/`` on sys.path). The locked Part 1 code is imported read-only:
``coessential`` (external/part1_locked/mvp/src) for CoessentialityModel and ``src.baselines`` (external/part1_locked)
for the GLS estimator, via :mod:`part2.part1`. Firewall F2: nothing in this package reads hu.MAP files or Part 1
test-evaluation results (tests/test_part2_firewall.py greps for it).
"""
