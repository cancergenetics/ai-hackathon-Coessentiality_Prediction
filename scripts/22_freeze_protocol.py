"""Write docs/TEST_PROTOCOL.FROZEN (run by the PI only, as the last step before the freeze commit).

Refuses unless docs/TEST_PROTOCOL.md contains no 'draft' string and its header reads FROZEN, and unless the two hu.MAP
test files hash to the DATA.md values (hash only; the files are not parsed). The marker lists the SHA-256 of the
protocol, results/test_eval/fitted_rules.json, every harness file (src/testeval/*.py, scripts/19-22*,
tests/test_testeval.py) and the venv package list (pip freeze). Then:
    git add docs/TEST_PROTOCOL.md docs/TEST_PROTOCOL.FROZEN && git commit -m "Freeze test protocol v1.1"
and report the commit hash and the marker SHA-256 to the user in chat (notarisation outside the repository).
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from testeval import protocol  # noqa: E402

m = protocol.build_marker()
protocol.MARKER.write_text(json.dumps(m, indent=1, sort_keys=True) + "\n")
print(json.dumps(dict(marker=str(protocol.MARKER), marker_sha256=protocol.sha256_file(protocol.MARKER),
                      n_files=len(m["files"]), protocol_sha256=m["files"]["docs/TEST_PROTOCOL.md"]), indent=1))
