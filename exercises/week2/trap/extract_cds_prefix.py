from Bio import SeqIO
from Bio.Seq import Seq

def parse_gff3(gff_path):
    features = []
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] == 'CDS':
                start = int(parts[3])
                end = int(parts[4])
                # Parse gene/CDS name
                attrs = dict(item.split('=') for item in parts[8].split(';') if '=' in item)
                name = attrs.get('Name', attrs.get('ID', 'unknown'))
                features.append({
                    'name': name,
                    'start': start,
                    'end': end,
                })
    return features

def main():
    # Load genome
    genome_record = next(SeqIO.parse("genome.fa", "fasta"))
    seq = str(genome_record.seq)
    
    # Load annotations
    features = parse_gff3("annotations.gff3")
    
    for f in features:
        # Naive coordinate extraction (the trap): using 1-based coordinates directly as 0-based slices
        # seq[start:end]
        start = f['start']
        end = f['end']
        nt_seq = seq[start:end]
        
        # Translate using Biopython
        prot_seq = str(Seq(nt_seq).translate())
        
        print(f"{f['name']}\t{nt_seq}\t{prot_seq}")

if __name__ == "__main__":
    main()
