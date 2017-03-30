

##########  test strand breaks composition : strand_breaks_composition() ##########


file <- "../data/Fu_2016/AfontovaGora3.archaic.q30.csv"
strand_breaks_composition(file, CtoT = TRUE, dist_from_ends = 20)

file <- "../data/moderns_lite/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.q30.csv"
strand_breaks_composition(file, CtoT = TRUE, dist_from_ends = 20)


