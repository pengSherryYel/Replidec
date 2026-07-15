#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_ROOT="${REPLIDEC_TEST_OUTPUT_DIR:-$(mktemp -d "${TMPDIR:-/tmp}/replidec-test.XXXXXX")}"
REPLIDEC_COMMAND="${REPLIDEC_COMMAND:-Replidec}"

mkdir -p "${OUTPUT_ROOT}"

for executable in "${REPLIDEC_COMMAND}" mmseqs hmmsearch blastp prodigal; do
    if ! command -v "${executable}" >/dev/null 2>&1; then
        echo "Required executable is unavailable: ${executable}" >&2
        exit 1
    fi
done

assert_summary() {
    local summary_file="$1"
    local minimum_rows="$2"

    if [[ ! -s "${summary_file}" ]]; then
        echo "Missing or empty summary: ${summary_file}" >&2
        exit 1
    fi

    if ! head -n 1 "${summary_file}" | grep -q $'^sample_name\t'; then
        echo "Unexpected summary header: ${summary_file}" >&2
        exit 1
    fi

    local rows
    rows="$(tail -n +2 "${summary_file}" | sed '/^[[:space:]]*$/d' | wc -l | tr -d ' ')"
    if (( rows < minimum_rows )); then
        echo "Expected at least ${minimum_rows} result rows in ${summary_file}, found ${rows}" >&2
        exit 1
    fi
}

"${REPLIDEC_COMMAND}" --version

cd "${SCRIPT_DIR}"

"${REPLIDEC_COMMAND}" \
    -p genome_table \
    -i example/genome_test.small.index \
    -w "${OUTPUT_ROOT}/genome_table" \
    -t 2
assert_summary "${OUTPUT_ROOT}/genome_table/prediction_summary.tsv" 4

"${REPLIDEC_COMMAND}" \
    -p multi_fasta \
    -i example/test.contig.small.fa \
    -w "${OUTPUT_ROOT}/multi_fasta" \
    -t 2
assert_summary "${OUTPUT_ROOT}/multi_fasta/prediction_summary.tsv" 1

"${REPLIDEC_COMMAND}" \
    -p protein_table \
    -i example/example.small.list \
    -w "${OUTPUT_ROOT}/protein_table" \
    -t 2
assert_summary "${OUTPUT_ROOT}/protein_table/prediction_summary.tsv" 4

echo "RepliDec integration tests passed. Outputs: ${OUTPUT_ROOT}"
