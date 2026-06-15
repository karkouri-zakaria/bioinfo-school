"""Small intentionally broken FASTA summarizer for Week 2.

Use an agentic IDE to fix this script one issue at a time. The intended
behavior is: read a FASTA file, print one row per sequence, and report sequence
length plus GC percentage.
"""

from pathlib import Path


def read_fasta(path):
    records = {}
    current_name = None
    current_seq = []

    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line.startswith(">"):
            if current_name:
                records[current_name] = "".join(current_seq)
            current_name = line[1:]
            current_seq = []
        else:
            current_seq.append(line)
    
    if current_name:
        records[current_name] = "".join(current_seq)

    return records


def gc_percent(sequence):
    gc = sequence.count("G") + sequence.count("C")
    return gc / len(sequence)


def main():
    records = read_fasta("example.fa")

    for name, sequence in records.items():
        print(name, len(sequence), gc_percent(sequence))


if __name__ == "__main__":
    main()
