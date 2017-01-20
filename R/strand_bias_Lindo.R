

#########  strand bias in Lindo data #######################

lindo_ancient <- get(load("../processed_data/lindo2016ancients-counts-table-strand-flank.rda"))
lindo_ancient_club <- club_signature_counts(lindo_ancient)

temp_data <- lindo_ancient_club[,-c(grep("left", colnames(lindo_ancient_club)), grep("right", colnames(lindo_ancient_club)))]

minus_indices <- grep("_-_", colnames(temp_data))
lindo_ancient_club_minus <- temp_data[,minus_indices]
lindo_ancient_club_plus <- temp_data[,-minus_indices]

filtered_data_minus <- filter_signatures_only_location(lindo_ancient_club_minus, max_pos = 20)
filtered_data_plus <- filter_signatures_only_location(lindo_ancient_club_plus, max_pos = 20)
filtered_data_plus <- filtered_data_plus[,-1]
filtered_data_minus <- filtered_data_minus[,-1]

sum(filtered_data_minus)/(sum(filtered_data_plus)+ sum(filtered_data_minus))


plot(log(filtered_data_minus[1,]), ylim=c(0,max(log(filtered_data_minus))), col="blue", type="l")
for(l in 2:25){
  lines(log(filtered_data_minus[l,]), col="blue")
}

for(l in 1:25){
  lines(log(filtered_data_plus[l,]), col="red")
}


sherpa_ancient <- get(load("../processed_data/sherpa2017-strand-flank.rda"))
sherpa_ancient_club <- club_signature_counts(sherpa_ancient)

temp_data <- sherpa_ancient_club[,-c(grep("left", colnames(sherpa_ancient_club)), grep("right", colnames(sherpa_ancient_club)))]

minus_indices <- grep("_-_", colnames(temp_data))
lindo_ancient_club_minus <- temp_data[,minus_indices]
lindo_ancient_club_plus <- temp_data[,-minus_indices]

filtered_data_minus <- filter_signatures_only_location(lindo_ancient_club_minus, max_pos = 20)
filtered_data_plus <- filter_signatures_only_location(lindo_ancient_club_plus, max_pos = 20)
filtered_data_plus <- filtered_data_plus[,-1]
filtered_data_minus <- filtered_data_minus[,-1]

sum(filtered_data_minus)/(sum(filtered_data_plus)+ sum(filtered_data_minus))

plot(log(filtered_data_minus[1,]), ylim=c(0,max(log(filtered_data_minus))), col="blue", type="l")
for(l in 2:5){
  lines(log(filtered_data_minus[l,]), col="blue")
}

for(l in 1:5){
  lines(log(filtered_data_plus[l,]), col="red")
}

mat <- lindo_ancient_club


filtered_data <- filter_signatures_only_location(lindo_ancient_club, max_pos = 20)

mat <- lindo_ancient_club
max_pos = 20

x <- as.character(colnames(mat))[1]
