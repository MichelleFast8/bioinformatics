import sys
from itertools import product

RED = "\033[31m"
GREEN = "\033[32m"
YELLOW = "\033[33m"
BLUE = "\033[34m"
RESET = "\033[0m"
COLORS = [RED, GREEN, YELLOW, BLUE]

def acid_to_codon(acid):
    codons = {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG"],
            "Y": ["TAT", "TAC"],
            "X": ["TAA", "TAG", "TGA"], # termination codon??
            "C": ["TGT", "TGC"],
            "W": ["TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "S": ["AGT", "AGC"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GGG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            }
    if acid not in codons.keys():
        print(f"Invalid amino acid abbreviation: {acid}")

    else:
        return codons[acid]

'''
Find all possible DNA sequences that could have sourced this polypeptide
'''
def main():	
    polypep_seq = list(sys.argv[1])

    # Retrieve and store each acid's associated codons, in order
    codon_lists = [acid_to_codon(acid) for acid in polypep_seq]

    # Calculate all possible DNA sequences
    dna_seqs = ["".join(codons) for codons in product(*codon_lists)]
    
    # Print Results
    print(f"Number of DNA Sequences: {len(dna_seqs)}")
    print("Polypeptide: ", end="")
    pretty_print(polypep_seq, False)
    for seq in dna_seqs:
        pretty_print(seq, True)
	

'''
Print codons in same color as corresponding acid
'''
def pretty_print(seq, is_dna):
    if (is_dna):
        split_seq = [seq[i:i+3] for i in range(0, len(seq), 3)]
        i=0
        for codon in split_seq:
            if i >= len(COLORS):
                i = 0
            print(f"{COLORS[i]}{codon}", end="")
            i = i+1
        print()
    else:
        i=0
        for acid in seq:
            if i >= len(COLORS):
                i = 0
            print(f"{COLORS[i]}{acid}", end="")
            i += 1
        print()


if __name__=="__main__":
	main()
