#!/usr/bin/env python3
# calc_metrics.py
# Alignment-block level performance metric calculation

import os
import sys
import re
from collections import defaultdict, Counter

# High threshold to accommodate large/megaplasmids assemblies
PLASMID_MAX_SIZE_THRESHOLD = 350_000

# --- Interval utilities ----------------------------------------------------

def merge_intervals(intervals):
    if not intervals:
        return 0
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return sum(e - s for s, e in merged)


# --- FASTA parsing ----------------------------------------------------------

def parse_fasta_regions(fasta_path):
    lengths          = {}
    is_plasmid       = {}
    is_circular      = {}
    circular_lengths = set()
    
    current_id     = None
    current_len    = 0
    current_header = ""

    if not fasta_path or not os.path.exists(fasta_path):
        return lengths, is_plasmid, is_circular, circular_lengths

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id is not None:
                    lengths[current_id] = current_len
                    hl = current_header.lower()
                    is_plasmid[current_id] = 'plasmid' in hl
                    
                    is_circ = bool(
                        re.search(r'circular\s*=\s*true', hl) or
                        'topology=circular' in hl or
                        'circular=1' in hl or
                        'circular=yes' in hl or
                        re.search(r'\bcircular\b', hl)
                    )
                    is_circular[current_id] = is_circ
                    if is_circ and current_len > 0:
                        circular_lengths.add(current_len)

                current_header = line[1:]
                current_id     = current_header.split()[0]
                current_len    = 0
            else:
                current_len += len(line)

    if current_id is not None:
        lengths[current_id] = current_len
        hl = current_header.lower()
        is_plasmid[current_id] = 'plasmid' in hl
        is_circ = bool(
            re.search(r'circular\s*=\s*true', hl) or
            'topology=circular' in hl or
            'circular=1' in hl or
            'circular=yes' in hl or
            re.search(r'\bcircular\b', hl)
        )
        is_circular[current_id] = is_circ
        if is_circ and current_len > 0:
            circular_lengths.add(current_len)

    return lengths, is_plasmid, is_circular, circular_lengths


# --- Predicted contig parsing ------------------------------------------------

def parse_predicted_contigs(pred_tsv, is_assembler=False):
    predicted = {}
    dropped_count = 0
    dropped_bases = 0

    with open(pred_tsv) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            contig_id = parts[0]
            try:
                size = int(parts[1])
            except ValueError:
                size = 0

            if is_assembler and size > PLASMID_MAX_SIZE_THRESHOLD:
                dropped_count += 1
                dropped_bases += size
                continue

            predicted[contig_id] = size

    if dropped_count > 0:
        print(f"  [Filter] Dropped {dropped_count} assembler contigs exceeding "
              f"{PLASMID_MAX_SIZE_THRESHOLD // 1000}kb "
              f"({dropped_bases:,} total bp of likely chromosomal leak).", file=sys.stderr)

    return predicted


# --- PAF parsing -----------------------------------------------------------

def parse_paf(paf_file, predicted_plasmid_contigs,
              plasmid_ref_ids, min_identity, is_assembler=False):
    contig_tp_intervals   = defaultdict(list)
    ref_covered_intervals = defaultdict(list)
    total_matches = 0
    total_aln_len = 0

    with open(paf_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 11:
                continue

            contig_id = parts[0]
            q_len     = int(parts[1])
            q_start   = int(parts[2])
            q_end     = int(parts[3])
            ref_id    = parts[5]
            r_start   = int(parts[7])
            r_end     = int(parts[8])
            matches   = int(parts[9])
            aln_len   = int(parts[10])

            if not is_assembler:
                # Check direct or substring match for contig IDs
                if contig_id not in predicted_plasmid_contigs:
                    matched = False
                    for pred_id in predicted_plasmid_contigs:
                        if contig_id == pred_id or contig_id.endswith(f"_{pred_id}") or pred_id.endswith(f"_{contig_id}"):
                            matched = True
                            break
                    if not matched:
                        continue
            else:
                if q_len > PLASMID_MAX_SIZE_THRESHOLD:
                    continue

            if matches == 0:
                continue

            block_identity = matches / aln_len if aln_len > 0 else 0.0
            if block_identity < (min_identity / 100.0):
                continue

            if ref_id not in plasmid_ref_ids:
                continue

            contig_tp_intervals[contig_id].append((q_start, q_end))
            ref_covered_intervals[ref_id].append((r_start, r_end))
            total_matches += matches
            total_aln_len += aln_len

    return (contig_tp_intervals, ref_covered_intervals,
            total_matches, total_aln_len)


# --- Circularization detection ---------------------------------------------------

def parse_circ_paf(circ_paf, pred_lengths, overlap_fraction=0.05):
    circular = defaultdict(list)
    try:
        with open(circ_paf) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 11:
                    continue
                qname = parts[0]
                tname = parts[5]

                if qname != tname:
                    continue

                seq_len = pred_lengths.get(qname, 0)
                if seq_len == 0:
                    continue

                min_overlap = seq_len * overlap_fraction
                q_start = int(parts[2])
                q_end   = int(parts[3])
                t_start = int(parts[7])
                t_end   = int(parts[8])

                end_to_start = (
                    q_end   >= seq_len - min_overlap and
                    t_start <= min_overlap
                )
                start_to_end = (
                    q_start <= min_overlap and
                    t_end   >= seq_len - min_overlap
                )

                if end_to_start or start_to_end:
                    circular[qname].append((min(q_start, t_start), max(q_end, t_end)))
    except FileNotFoundError:
        pass

    return circular


# --- Size category ------------------------------------------------------------

def get_size_category(size):
    if size < 1_000:       return '0-1kb'
    elif size < 5_000:     return '1-5kb'
    elif size < 10_000:    return '5-10kb'
    elif size < 100_000:   return '10-100kb'
    else:                  return '>100kb'


# --- Main --------------------------------------------------------------------

def main():
    if len(sys.argv) not in (7, 8):
        print(
            "Usage: calc_metrics.py "
            "<paf_file> <ref_fasta> <predicted_contigs_tsv> "
            "<min_identity> <circ_paf> <plasmid_fasta> [fallback_size_cat]",
            file=sys.stderr
        )
        sys.exit(1)

    paf_file          = sys.argv[1]
    ref_fasta         = sys.argv[2]
    pred_tsv          = sys.argv[3]
    min_identity      = float(sys.argv[4])
    circ_paf          = sys.argv[5]
    plasmid_fasta     = sys.argv[6]
    fallback_size_cat = sys.argv[7] if len(sys.argv) > 7 else 'NA'

    missing_input = False
    if not os.path.exists(pred_tsv) or not os.path.exists(paf_file):
        missing_input = True

    # --- Step 1: Parse reference FASTA ------------------------------------------------
    ref_lengths, ref_is_plasmid, _, _ = parse_fasta_regions(ref_fasta)
    plasmid_ref_ids = {r for r, p in ref_is_plasmid.items() if p}

    # Fallback: If no header explicitly contained 'plasmid', treat ALL ref sequences as plasmid refs
    if not plasmid_ref_ids:
        plasmid_ref_ids = set(ref_lengths.keys())

    total_ref_size = sum(ref_lengths[r] for r in plasmid_ref_ids)

    # --- Step 2: Broad Assembler Detection ------------------------------------------
    # Check tool keywords or file paths for assembler indications
    combined_paths = (pred_tsv + " " + paf_file + " " + plasmid_fasta).lower()
    is_assembler = any(k in combined_paths for k in ["plassembler", "assembly", "assembler", "unpolished", "flye", "unicycle", "spades"])

    if not missing_input:
        predicted_plasmid_contigs = parse_predicted_contigs(pred_tsv, is_assembler=is_assembler)
    else:
        predicted_plasmid_contigs = {}

    total_predicted_size = sum(predicted_plasmid_contigs.values())

    # --- Step 3: Parse PAF at alignment-block level ---------------------------------
    if not missing_input and os.path.exists(paf_file):
        (contig_tp_intervals,
         ref_covered_intervals,
         total_matches,
         total_aln_len) = parse_paf(
            paf_file,
            predicted_plasmid_contigs,
            plasmid_ref_ids,
            min_identity,
            is_assembler=is_assembler
        )
    else:
        contig_tp_intervals = defaultdict(list)
        ref_covered_intervals = defaultdict(list)
        total_matches = 0
        total_aln_len = 0

    # --- Step 4: Base-pair level TP, FP, FN ------------------------------------------
    total_ref_covered = sum(merge_intervals(intervals) for intervals in ref_covered_intervals.values())

    tp_bases = total_ref_covered
    fn_bases = max(total_ref_size - total_ref_covered, 0)
    fp_bases = max(total_predicted_size - tp_bases, 0)

    # --- Step 5: Benchmarking metrics ------------------------------------------------
    bp_precision = (tp_bases / (tp_bases + fp_bases)
                    if (tp_bases + fp_bases) > 0 else 0.0)
    bp_recall    = (total_ref_covered / total_ref_size
                    if total_ref_size > 0 else 0.0)

    f1 = (2 * bp_precision * bp_recall / (bp_precision + bp_recall)
          if (bp_precision + bp_recall) > 0 else 0.0)

    completeness = round(bp_recall * 100, 4)
    contamination = round(
        fp_bases / total_predicted_size * 100
        if total_predicted_size > 0 else 0.0,
        4
    )
    identity = round(
        total_matches / total_aln_len * 100
        if total_aln_len > 0 else 0.0,
        4
    )
    accuracy = round(
        tp_bases / (tp_bases + fp_bases + fn_bases)
        if (tp_bases + fp_bases + fn_bases) > 0 else 0.0,
        4
    )

    bp_precision = round(bp_precision, 4)
    bp_recall    = round(bp_recall,    4)
    f1           = round(f1,           4)

    # --- Step 6: Contig-level TP/FP (50% rule) -----------------------------------
    tp_count = 0
    fp_count = 0

    if not is_assembler:
        for contig_id, size in predicted_plasmid_contigs.items():
            mapped_tp_bases = merge_intervals(contig_tp_intervals[contig_id])
            if size > 0 and (mapped_tp_bases / size) > 0.50:
                tp_count += 1
            else:
                fp_count += 1

    contig_precision = round(
        tp_count / (tp_count + fp_count)
        if (tp_count + fp_count) > 0 else 0.0,
        4
    )

    # --- Step 7: Detect Circular Contigs ---------------------------------------------
    ref_circular_intervals = defaultdict(list)

    if is_assembler:
        # For assemblers (e.g. Plassembler), all valid PAF alignments represent predicted plasmids.
        ref_circular_intervals = ref_covered_intervals
    else:
        circular_pred_intervals = parse_circ_paf(circ_paf, predicted_plasmid_contigs)
        circular_contig_ids = set(circular_pred_intervals.keys())
        header_circular_lengths = set()

        if plasmid_fasta and os.path.exists(plasmid_fasta):
            _, _, header_circular_map, header_circular_lengths = parse_fasta_regions(plasmid_fasta)
            for cid, is_circ in header_circular_map.items():
                if is_circ:
                    circular_contig_ids.add(cid)

        circular_sizes = set(header_circular_lengths)
        for cid in circular_contig_ids:
            if cid in predicted_plasmid_contigs:
                circular_sizes.add(predicted_plasmid_contigs[cid])

        if os.path.exists(paf_file):
            with open(paf_file) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 11:
                        continue
                    contig_id = parts[0]
                    q_len     = int(parts[1])
                    ref_id    = parts[5]
                    r_start   = int(parts[7])
                    r_end     = int(parts[8])
                    matches   = int(parts[9])
                    aln_len   = int(parts[10])

                    if contig_id not in circular_contig_ids and q_len not in circular_sizes:
                        continue

                    if ref_id not in plasmid_ref_ids or matches == 0:
                        continue

                    block_identity = matches / aln_len if aln_len > 0 else 0.0
                    if block_identity < (min_identity / 100.0):
                        continue

                    ref_circular_intervals[ref_id].append((r_start, r_end))

    circular_tp_bases = sum(merge_intervals(intervals) for intervals in ref_circular_intervals.values())

    circularization_pct = round(
        circular_tp_bases / total_ref_size * 100
        if total_ref_size > 0 else 0.0,
        4
    )

    # --- Step 8: Size category ------------------------------------------------
    if plasmid_ref_ids:
        cats     = [get_size_category(ref_lengths[r]) for r in plasmid_ref_ids]
        size_cat = Counter(cats).most_common(1)[0][0]
    else:
        size_cat = fallback_size_cat

    # --- Step 9: Output --------------------------------------------------------
    print(
        f"{tp_bases}\t{fp_bases}\t{fn_bases}\t{accuracy}\t"
        f"{circular_tp_bases}\t{circularization_pct}\t"
        f"{bp_recall}\t{contig_precision}\t{bp_precision}\t{f1}\t"
        f"{identity}\t{completeness}\t"
        f"{contamination}\t{size_cat}"
    )


if __name__ == '__main__':
    main()
