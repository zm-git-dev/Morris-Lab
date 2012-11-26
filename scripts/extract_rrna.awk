#!/bin/awk
#
# Extract a fragment of a FASTA sequence centered at a particular position.
# NOTE: this script requires a fasta file on stdin that 
# contains the named sequence to be extracted.
# 
# usage: 
#  extract_rrna.awk <sequence_name> <position_in_seuqnce> [<length_of_oligo>] 
#  e.g.: extract_rrna.awk 45S 12377 30 <rrna.fa
#
# This is used to extract oligos for experimental rRNA
# depletion. Typical usage is to first identify candidate regions by
# sorting the robosomal positions by read depth.  Then take the
# position with the greatest read depth and construct an oligo from
# the rRNA reference around that position.
#
# samtools view -S -b rrna_aligned.sam >rrna.bam
# samtools sort rrna.bam rrna_s
# samtools index rrna_s.bam
# samtools depth rrna_s.bam | sort -k 3rn >rrna_depth.txt
# head rrna_depth.txt
# <grab the highest depth position>
# awk -f extract_rrna.awk 45S 12733 <~morrislab/genome/mm9/rrna.fa
# --> TGTGGAACCTGGCGCTAAACCATTCGTAGA
#
function process_sequence(seqname, sequence, target, position) {
    n = split(seqname, names, "|");
    for (i = 1; i <= n; i++) {
	where = match(names[i], target ".*");
	if (where != 0) {
	    print substr(sequence, position-15, 30)
	    exit;
	}	    
    }	    
}

BEGIN {
    target = ARGV[1];
    position = ARGV[2] + 1;
    in_sequence = 0;
    delete ARGV[1]
    delete ARGV[2]
}

# skip comments in FASTA files.
#
/\w*\#/ { next; }

# beginning a sequence with '>name' line foloowed by nucleotide sequences
# split across multiple lines.
/^>/  {
    if (in_sequence == 1) {
	process_sequence(seqname, sequence, target, position);
	sequence = "";
    }


    seqname = substr($_, 2)
    in_sequence=1;
    next;
}

in_sequence == 1 {
    sequence = sequence $0;
}

END {
    if (in_sequence == 1) {
	process_sequence(seqname, sequence, target, position);
    }



}
