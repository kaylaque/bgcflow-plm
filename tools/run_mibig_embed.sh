#!/usr/bin/env bash
# Resiliently build resources/mibig_esm2/index.faiss on a shared machine.
#
# Why this exists: the embedding process has an unbounded CPU-memory-growth
# issue (RSS climbed to ~28GB after only ~500 of 39,992 proteins in a real
# run — see docs/systems/plm-novelty.md). Left unguarded on a shared box,
# that risks the kernel OOM killer picking victims across OTHER users'
# processes, not just ours. This wrapper caps our own job to a memory
# budget via a user-level systemd scope (no sudo needed) so if it grows too
# large, only OUR job gets killed cleanly — then retries, resuming from the
# last checkpoint shard (embed_sequences() in embed_esm2.py checkpoints by
# protein id, so a killed run never starts over from zero).
#
# Usage:
#   tools/run_mibig_embed.sh
#   MEMORY_MAX=24G MAX_RETRIES=100 tools/run_mibig_embed.sh
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."

MEMORY_MAX="${MEMORY_MAX:-16G}"
MAX_RETRIES="${MAX_RETRIES:-100}"
RETRY_DELAY_SECONDS="${RETRY_DELAY_SECONDS:-10}"
TARGET="resources/mibig_esm2/index.faiss"

echo "[run_mibig_embed] memory cap: ${MEMORY_MAX}, max retries: ${MAX_RETRIES}"

for attempt in $(seq 1 "${MAX_RETRIES}"); do
    if [[ -f "${TARGET}" ]]; then
        echo "[run_mibig_embed] ${TARGET} already exists — done."
        exit 0
    fi

    echo "[run_mibig_embed] attempt ${attempt}/${MAX_RETRIES} starting at $(date -Iseconds)"
    if systemd-run --user --scope --quiet \
        -p "MemoryMax=${MEMORY_MAX}" -p "MemorySwapMax=0" -- \
        snakemake --use-conda --cores 8 --rerun-triggers mtime -- "${TARGET}"; then
        echo "[run_mibig_embed] completed successfully."
        exit 0
    fi

    echo "[run_mibig_embed] attempt ${attempt} exited non-zero (likely the memory cap or an interruption) — retrying in ${RETRY_DELAY_SECONDS}s"
    sleep "${RETRY_DELAY_SECONDS}"
done

echo "[run_mibig_embed] exhausted ${MAX_RETRIES} attempts without producing ${TARGET}." >&2
exit 1
