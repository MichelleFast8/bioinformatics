import sys

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

def main():	
    polypep_seq = list(sys.argv[1])

    # Retrieve and store each acid's associated codons
    for i in range(len(polypep_seq)):
        polypep_seq[i] = (polypep_seq[i], acid_to_codon(polypep_seq[i]))

    print(polypep_seq)
	

if __name__=="__main__":
	main()
