"""Regression guardrail: plm_novelty.smk must never write into antiSMASH output dirs."""
import re
from pathlib import Path

SMK_PATH = Path(__file__).parents[2] / "workflow" / "rules" / "plm_novelty.smk"
FORBIDDEN_PATTERNS = [
    r"data/interim/antismash",
    r"data/processed/\{name\}/antismash",
    r"data/processed/\{name\}/(tables|bigslice|bigscape|roary|automlst|gecco)",
]


def test_smk_file_exists():
    assert SMK_PATH.exists(), f"plm_novelty.smk not found at {SMK_PATH}"


def test_no_antismash_output_writes():
    """
    Verify that plm_novelty.smk contains no 'output:' block that writes
    into the antiSMASH interim or processed directories.
    """
    content = SMK_PATH.read_text()

    # Parse output blocks: find 'output:' sections and check their paths
    in_output_block = False
    output_lines = []
    for line in content.splitlines():
        stripped = line.strip()
        if stripped.startswith("output:"):
            in_output_block = True
            continue
        if in_output_block:
            # End of output block: next non-indented keyword
            if re.match(r"^\s*(input|params|conda|log|shell|run|threads|resources|wildcard_constraints)\s*:", line):
                in_output_block = False
                continue
            output_lines.append(line)

    output_text = "\n".join(output_lines)

    for pattern in FORBIDDEN_PATTERNS:
        matches = re.findall(pattern, output_text)
        assert not matches, (
            f"plm_novelty.smk output block contains forbidden path pattern '{pattern}'. "
            f"The PLM lane must not write into antiSMASH directories."
        )
