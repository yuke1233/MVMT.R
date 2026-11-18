#!/usr/bin/env bash
set -euo pipefail

BED="../chr22_gene.bed"

if [[ ! -f "$BED" ]]; then
    echo "ERROR: Cannot find $BED"
    exit 1
fi

# Read target genes from BED column 4
declare -A keep

while IFS= read -r gene; do
    [[ -n "$gene" ]] && keep["$gene"]=1
done < <(awk 'NF >= 4 && $1 !~ /^#/ {print $4}' "$BED" | sort -u)

echo "Target genes: ${#keep[@]}"

# Preview by default; pass --delete to actually remove files
delete_mode=0
if [[ "${1:-}" == "--delete" ]]; then
    delete_mode=1
fi

removed=0
retained=0

for f in *.ldm.full.bin *.ldm.full.info; do
    [[ -f "$f" ]] || continue

    gene="${f%.ldm.full.bin}"
    gene="${gene%.ldm.full.info}"

    if [[ -n "${keep[$gene]+yes}" ]]; then
        retained=$((retained + 1))
    else
        if (( delete_mode )); then
            rm -- "$f"
            echo "DELETED: $f"
        else
            echo "TO DELETE: $f"
        fi
        removed=$((removed + 1))
    fi
done

echo
echo "========== Summary =========="
echo "Target genes: ${#keep[@]}"
echo "Files to retain: $retained"
echo "Files to delete: $removed"
if (( ! delete_mode )); then
    echo "Preview only. Run with --delete to apply."
fi
echo "============================="
