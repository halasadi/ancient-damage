

##############  test flanking base composition  ########################

data <- read.csv("../data/Fu_2016/AfontovaGora3.archaic.q30.csv", header=FALSE)
flanking_base_composition(file="../data/Fu_2016/AfontovaGora3.archaic.q30.csv",
                          pattern = "C->A")
