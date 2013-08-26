source("../morrislib.R")
source("profiling.R")

#gene="NR_029560"   # Mir150
#gene='NM_TEST'		# 'test'
gene='NM_144903'	# Aldob '-'
gene='NM_011434'	# 'Sod1'
gene='NM_013541'	# 'Gstp1'
gene='NM_001005419'	# 'Ado'  single exon '-'
gene="NM_007409"   # ADH1
gene="NR_029600"   # Mir122a   single exon '+'
gene="NM_133862"    ## Fgg
gene='NM_016978'	# Oat '-'
gene='NM_009790'    ## CALM1 - Calmodulin - no B-sheets
gene ="NM_009654" ## Alb
gene='NM_011044'	# 'Pck1'
gene='NM_016978'    # Oat '-'
gene='NM_001005419'    # 'Ado'  single exon '-'
gene="NM_007409"   # ADH1
gene="RPS23"  

plot.new()

kg <- morris.getknowngenes(genome, gene=gene, group=NULL)
rownames(kg) <- kg$name
stopifnot(nrow(kg) == 1)
        
refseq = kg[1,'name']


df = morris.getalignments("061113_A", refseq)

p = profile(df, kg[refseq,])

print(p)
rpos = plot(p, minlen=0, units="aa")

