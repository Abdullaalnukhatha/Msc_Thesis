#!/usr/bin/env python3
# get_ref_size_category.py
# Calculate size category from FASTA file

import sys

def get_size_category(size):
    if size < 1_000:
        return '0-1kb'
    elif size < 5_000:
        return '1-5kb'
    elif size < 10_000:
        return '5-10kb'
    elif size < 100_000:
        return '10-100kb'
    else:
        return '>100kb'

def main():
    fasta = sys.argv[1]

    sizes = []
    current_len = 0
    with open(fasta) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_len > 0:
                    sizes.append(current_len)
                current_len = 0
            else:
                current_len += len(line)
        if current_len > 0:
            sizes.append(current_len)

    if not sizes:
        print('NA')
        return

    # Dominant size category (most common among plasmids in this isolate)
    from collections import Counter
    cats = [get_size_category(s) for s in sizes]
    dominant = Counter(cats).most_common(1)[0][0]
    print(dominant)

if __name__ == '__main__':
    main()
