# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "click",
#     "numpy",
# ]
# ///
import random
import numpy as np
import click

seed = 12345

random.seed(seed)
np.random.seed(seed)

sequences_to_avoid = ["gtgtatcggatgtcagttgc",
                      "gtataatgcagacctgctgcCGGATGAAAACGAGA",
                      "AATGATACGGCGACCACCGAGATCTACAC",
                      "ACACTCTTTCCCTACACGACGCTCTTCCGATCTGTGTATCGGATGTCAGTTGC",
                      "CAAGCAGAAGACGGCATACGAGAT",
                      "GTGACTGGAGTTCCTTGGCACCCGAGAATTCCATCTCGTTTTCATCCGGCAGC",
                      "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGAgatgatgattacTAATACGACTCACTATAACTGGAAG",
                      "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACActcgttttcatccggcag"]

# Constraints
def gc_content(seq: str) -> float:
    """Return GC percentage."""
    gc = seq.count("G") + seq.count("C")
    return 100 * gc / len(seq)

def check_uniqueness(seq: str, sequences_to_avoid: list) -> bool:
    """Check whether sequence is in list of sequences to avoid"""
    unique = True
    for seq_avoid in sequences_to_avoid:
        if seq.lower() in seq_avoid.lower() or seq.lower() in reverse_complement(seq_avoid).lower():
            unique = False
    return unique

def hamming_distance(a: str, b: str) -> int:
    """Compute Hamming distance between sequences"""
    if len(a) == len(b):
        a = list(a)
        b = list(b)
        mismatches = [x[0] != x[1] for x in zip(a, b)]
        hamming_dist = sum(mismatches)
        return hamming_dist
    else:
        return

def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    seq = seq.upper()
    complement = str.maketrans("ATGC", "TACG")
    return seq.translate(complement)[::-1]

def write_barcodes_as_fasta(barcodes: list, outfile: str, prefix="") -> str:
    """Writes a list of barcodes out to a file in fasta format"""
    if prefix == "":
        prefix = "BARCODE"
    with open(outfile, "w") as out:
        for i, barcode in enumerate(barcodes):
            out.write(f">{prefix}{i+1:04d}\n{barcode}\n")
    return outfile

def generate_barcodes(number,
                      length,
                      hamming,
                      min_gc,
                      max_gc):
    barcodes = []
    gc = random.randint(min_gc,max_gc)
    while len(barcodes) < number:
        number_GCs =  int((gc / 100) * length)
        number_ATs = length - number_GCs
        prop_Gs = float(np.random.normal(0.5, 0.1, 1)[0])
        number_Gs = int(prop_Gs * number_GCs)
        number_Cs = number_GCs - number_Gs
        prop_As = float(np.random.normal(0.5, 0.1, 1)[0])
        number_As = int(number_ATs * prop_As)
        number_Ts = number_ATs - number_As
        barcode = "G" * number_Gs + "C" * number_Cs + "A" * number_As + "T" * number_Ts
        barcode = list(barcode)
        random.shuffle(barcode)
        barcode = "".join(barcode)
        if reverse_complement(barcode) != barcode:
            if all([hamming_distance(barcode, b) >= hamming for b in barcodes]):
                rev_comp_barcodes = list([reverse_complement(b) for b in barcodes])
                if all([hamming_distance(barcode, b) >= hamming for b in rev_comp_barcodes]):
                    if check_uniqueness(barcode, sequences_to_avoid):
                        barcodes.append(barcode)
    return barcodes

@click.command()
@click.option("-n", "--number", type=click.INT, default=100)
@click.option("-l", "--length", type=click.INT, default=16)
@click.option("-h", "--hamming", type=click.INT, default=6)
@click.option("-m", "--min_gc", type=click.INT, default=40)
@click.option("-a", "--max_gc", type=click.INT, default=60)
@click.option("-o", "--outfile", type=click.STRING, default="barcodes.fasta")
def cli(number, length, hamming, min_gc, max_gc, outfile):
    barcodes = generate_barcodes(number, length, hamming, min_gc, max_gc)
    filename = write_barcodes_as_fasta(barcodes, outfile)
    print(f"{len(barcodes)} barcodes written to {filename}")

if __name__ == "__main__":
    cli()