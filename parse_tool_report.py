#!/usr/bin/env python3
# parse_tool_report.py
# Parse tool-specific contig classification report using pandas

import sys
import os
import pandas as pd
from pathlib import Path


def get_fasta_lengths(fasta_path):
    lengths = {}
    current_id = None
    current_len = 0
    try:
        with open(fasta_path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_id is not None:
                        lengths[current_id] = current_len
                    current_id = line[1:].split()[0]
                    current_len = 0
                else:
                    current_len += len(line)
        if current_id is not None:
            lengths[current_id] = current_len
    except FileNotFoundError:
        pass
    return lengths


def parse_mobsuite(report_file, plasmid_fasta=None):
    if os.path.getsize(report_file) == 0:
        return []
    try:
        df = pd.read_csv(report_file, sep='\t')
        if df.empty or 'molecule_type' not in df.columns:
            return []
        df_filtered = df[df['molecule_type'].str.strip().str.lower() == 'plasmid']
        return list(zip(df_filtered['contig_id'].str.strip(), df_filtered['size'].fillna(0).astype(int)))
    except Exception:
        return []


def parse_plasmer(report_file, plasmid_fasta=None):
    if os.path.getsize(report_file) == 0:
        return []
    if plasmid_fasta is not None:
        SRR = os.path.basename(report_file).split(".")[0].strip()
        plasmid_fasta = os.path.join(os.path.dirname(report_file), f"{SRR}.plasmer.predPlasmids.fa")
    fasta_lengths = get_fasta_lengths(plasmid_fasta) if plasmid_fasta else {}

    try:
        df = pd.read_csv(report_file, sep='\t', header=None, names=['contig_id', 'molecule_type'])
        if df.empty:
            return []
        df_filtered = df[df['molecule_type'].str.strip().str.lower() == 'plasmid']
        
        plasmids = []
        for cid in df_filtered['contig_id'].str.strip():
            plasmids.append((cid, fasta_lengths.get(cid, 0)))
        return plasmids
    except Exception:
        return []


def parse_plassembler(report_file, plasmid_fasta=None):
    if os.path.getsize(report_file) == 0:
        return []
    try:
        df = pd.read_csv(report_file, sep='\t')
        if df.empty or 'contig' not in df.columns:
            return []
        df = df.dropna(subset=['contig'])
        return list(zip(df['contig'].str.strip(), df['length'].fillna(0).astype(int)))
    except Exception:
        return []


def parse_plasmidec(report_file, plasmid_fasta=None):
    if os.path.getsize(report_file) == 0:
        return []
    if plasmid_fasta is not None:
        plasmid_fasta = os.path.join(os.path.dirname(report_file), "plasmid_contigs.fasta")
    fasta_lengths = get_fasta_lengths(plasmid_fasta) if plasmid_fasta else {}

    try:
        df = pd.read_csv(report_file)
        if df.empty or 'Combined_prediction' not in df.columns:
            return []
        df['pred_str'] = df['Combined_prediction'].astype(str).str.strip().str.lower()
        df_filtered = df[df['pred_str'].isin(['1', '1.0', 'plasmid'])]

        plasmids = []
        for cid in df_filtered['Contig_name'].str.strip():
            plasmids.append((cid, fasta_lengths.get(cid, 0)))
        return plasmids
    except Exception:
        return []


def parse_platon(report_file, plasmid_fasta=None):
    if os.path.getsize(report_file) == 0:
        return []
    try:
        df = pd.read_csv(report_file, sep='\t')
        if df.empty or 'ID' not in df.columns:
            return []
        df = df.dropna(subset=['ID'])
        return list(zip(df['ID'].str.strip(), df['Length'].fillna(0).astype(int)))
    except Exception:
        return []


def parse_rfplasmid(report_file, plasmid_fasta=None):
    if os.path.getsize(report_file) == 0:
        return []
    try:
        # Read file. We let pandas infer separator (handles commas or spaces smoothly)
        # and parse the first column dynamically since R row names often lack a header name.
        df = pd.read_csv(report_file, sep=None, engine='python')
        
        # Filter for rows predicted as plasmids ('p')
        df_filtered = df[df['prediction'].astype(str).str.strip().str.lower() == 'p']
        
        # Return a list of tuples: (contigID, contig_length)
        return list(zip(
            df_filtered['contigID'].astype(str).str.strip(), 
            df_filtered['contig_length'].fillna(0).astype(int)
        ))
    except Exception as e:
        print(f"RFPlasmid parse error: {e}", file=sys.stderr)
        return []


PARSERS = {
    'MOB-suite':  parse_mobsuite,
    'Plasmer':    parse_plasmer,
    'Plassembler': parse_plassembler,
    'PlasmidEC':  parse_plasmidec,
    'Platon':     parse_platon,
    'RFPlasmid':  parse_rfplasmid,
}


def main():
    if len(sys.argv) < 3:
        print("Usage: parse_tool_report.py <tool_name> <report_file> [plasmid_fasta]",
              file=sys.stderr)
        sys.exit(1)

    tool = sys.argv[1]
    report_file = sys.argv[2]
    plasmid_fasta = sys.argv[3] if len(sys.argv) > 3 else None

    if tool not in PARSERS:
        print(f"ERROR: Unknown tool '{tool}'. "
              f"Valid options: {list(PARSERS.keys())}", file=sys.stderr)
        sys.exit(1)

    if not os.path.isfile(report_file):
        print(f"ERROR: Report file not found: {report_file}", file=sys.stderr)
        sys.exit(1)

    plasmids = PARSERS[tool](report_file, plasmid_fasta)

    for contig_id, size in plasmids:
        print(f"{contig_id}\t{size}")


if __name__ == '__main__':
    main()