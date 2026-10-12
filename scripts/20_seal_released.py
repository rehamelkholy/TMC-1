"""Released-regime evaluation (primary, sensitivities A4(a)/A4(b), v3 feedback learner).

Thin command-line wrapper added for the release: the evaluation itself lives in src/seal_released.py.

Usage: python scripts/20_seal_released.py <output.json>
Reads data/frozen/{candidates,lee_list_entry_labels_purged,post2021_detections}.csv; writes the JSON and a
.sha256 sidecar. The sealed reference copy is results/sealed/released.json.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from seal_released import main  # noqa: E402

if __name__ == "__main__":
    main()
