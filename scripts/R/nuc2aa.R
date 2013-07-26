require(Biostrings)
require(genomes)



Sys.setenv(email='cswarth@gmail.com')


nuc2aa <- function(gene) {

    ## use bioconductor - ShortRead package to read the fasta format.
    ## http://stackoverflow.com/a/9331406/1135316
    ## http://www.bioconductor.org/packages/release/bioc/html/ShortRead.html

    ## looks like BioStrings is a better package to read fasta strings.
    ## use the genomes package the fetch data from NCBI
    ## http://www.bioconductor.org/packages//2.10/bioc/manuals/genomes/man/genomes.pdf
    ## Results must be saved to a temp file
    ## http://www.biostars.org/p/75700/#75709
    tmp = tempfile()
    efetch(gene, showURL=FALSE, db="nucleotide", retmode="text", rettype="fasta", destfile=tmp, complexity=3, strand=1)
    records = readBStringSet(tmp)
    mrna = as(records[[1]], "DNAString")
    protein = as(records[[2]], "AAString")

    ## http://a-little-book-of-r-for-bioinformatics.readthedocs.org/en/latest/src/chapter7.html
    ## Find all occurrences of codon "codon" in sequence "sequence"
    match.atg <- matchPattern("atg", mrna)

    ## Two equivalent ways to extract a view as an XString object:
    s2b <- subseq(subject(match.atg), start=start(match.atg)[1], end=length(mrna))
    codons = as.character(codons(s2b))
    aa = strsplit(as.character(translate(s2b)), '')[[1]]
    n = match('*', aa)

    data.frame(codon=codons[1:n], aa=aa[1:n])

}

nuc2aa("NM_009790")

