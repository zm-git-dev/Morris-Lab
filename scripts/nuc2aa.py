"""
    Print codon and amino acid sequences for genes.

    Dave needs tables of the form "ATG   M" for genes under study.
    This script will take the refeseq name of a
    gene and download the mRNA nucleotide sequence from NCBI.  This
    script calculates its own amino acid sequence based on the
    nucleotide sequece, but it also verifies that this matches the
    peptide supplied by NCBI.


"""


import sys
import getopt
import collections

from Bio import Entrez, SeqIO


import sys
import getopt
import collections


# A nice fast way to initialize a codon table in python.
# Borrowed from PW Collingridge
# http://www.petercollingridge.co.uk/python-bioinformatics-tools/codon-table
bases = ['T', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))

def nt2protein(s1):
    """
    Return a protein sequence for a given DNA sequence.
    """
    # Look for start codon
    # Read off proteins until stop codon.
    n = s1.find("ATG")
    peptide = []
    for i in range(n,len(s1),3):
        codon = s1[i:i+3]
        if codon in codon_table:
            aa = codon_table[codon]
            if aa == "*":
                break
            peptide.append((codon,aa))
        else:
            print "Found non-canonical codon \"%s\"" % (codon)

    return (peptide)


def main(argv = None):
    Entrez.email = "warthc@uw.washington.edu"     # Always tell NCBI who you are

    genes = (
        'NM_009654', 
        'NM_009692', ## Apoa1 apolipoprotein A-I
        'NM_009693', ## Mus musculus apolipoprotein B (Apob)
        'NM_181849', ## FGB
        'NM_133862'  ## Fgg
      )
    for gene in genes:
        print gene
        handle = Entrez.efetch(db="nucleotide", 
                               id=gene,
                               complexity=3,
                               rettype="fasta", 
                               strand=1)
        records = list(SeqIO.parse(handle, "fasta"))
        mrna = records[0]
        protein = records[1]
        codons = nt2protein(str(mrna.seq))
        assert("".join([p for (c,p) in codons]) == str(protein.seq))
        handle.close()

        n1 = mrna.description.find('(')
        n2 = mrna.description.find(')')
        name = mrna.description[n1+1:n2]
        with open(name+".csv", "w") as f:
            for (c,p) in codons:
                f.write("{0}\t{1}\n".format(c,p))
            f.close()

        
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

