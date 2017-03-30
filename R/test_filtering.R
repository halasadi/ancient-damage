

mat <- get(load("../data/Fu_2016/Fu_2016.rda"))
max_pos = 20
input_pos <- 1:max_pos;

flanking_bases <- 1
x <- as.character(colnames(mat))[1]
mutations <- as.character(sapply(as.character(colnames(mat)), function(x) return (paste0(strsplit(x, "")[[1]][(flanking_bases + 1): (flanking_bases + 4)], collapse=""))))
pos <- as.numeric(sapply(as.character(colnames(mat)), function(x) return (strsplit(x,"_")[[1]][4])))

mutations_pos <- paste0(mutations, "_", pos)

mat_reduced <- as.numeric()
for(l in 1:dim(mat)[1]){
  mat_reduced <- rbind(mat_reduced, tapply(mat[l,], mutations_pos,  sum))
}

rownames(mat_reduced) <- rownames(mat)

####################  filter by mutation flank ############################

x <- as.character(colnames(mat))[1]
mat_reduced <- filter_by_mutation_flank(mat)

####################   filter by mutations ############################

mat_reduced <- filter_by_mutation(mat)

###################   filter by C->T and pos #############################

mat <- get(load("../data/Fu_2016/Fu_2016.rda"))
max_pos = 20
mat_reduced <- filter_by_pos_pattern(mat)
mat_reduced <- filter_by_pos_pattern(mat, pattern = "C->A")

###################  filter by mutation + flanking base + position ###########

mat <- get(load("../data/Fu_2016/Fu_2016.rda"))
max_pos = 20
mat_reduced <- filter_by_mutation_flank_pos(mat)

################  filter out strand ########################

mat <- get(load("../data/Fu_2016/Fu_2016.rda"))
max_pos = 20
mat_reduced <- filter_out_strand(mat)

#############  filter out strand break  ######################

mat <- get(load("../data/Fu_2016/Fu_2016.rda"))
max_pos = 20
mat_reduced <- filter_out_strand_break(mat)

############  filter out by pos   #########################

mat <- get(load("../data/Fu_2016/Fu_2016.rda"))
max_pos = 20
mat_reduced <- filter_by_pos(mat)



