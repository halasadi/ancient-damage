names = c("AfontovaGora3", "BerryAuBac", "Bockstein", "Chaudardes1",
          "Cioclovina1", "Continenza", "ElMiron", "GoyetQ116-1",
          "GoyetQ376-19", "GoyetQ53-1", "GoyetQ56-16", "HohleFels49",
          "Kostenki12", "Kostenki14", "KremsWA3", "Muierii",
          "Muierii2", "Ofnet", "Ostuni1", "Ostuni2",
          "Paglicci108", "Paglicci133", "Pavlov1", "Ranchot88",
          "Rigney1", "Rochedane", "Vestonice13", "Vestonice14",
          "Vestonice15", "Vestonice16", "Vestonice43", "Villabruna")

UDG = c("noUDG", "noUDG", "noUDG", "noUDG", 
        "noUDG", "UDG", "UDG", NA, 
        "noUDG", "noUDG", "noUDG", "noUDG", 
        "UDG", "UDG", NA, NA, 
        "UDG", "noUDG", "UDG", "UDG", 
        "noUDG", "noUDG", "UDG", "noUDG", 
        "noUDG", NA, "UDG", "UDG", 
        "UDG", "UDG", "UDG", "noUDG")
UDG = as.factor(UDG)
MT = c(0.99, 0.97,0.97, 0.93, 
       0.88, 0.99, 0.99, 0.95,
       0.92, 0.82, 0.84, 0.99, 
       0.99, 0.99, NA, NA, 
       0.93, 1, 0.82, 0.83, 
       0.93, 0.83, 0.99, 0.99, 
       0.9, 0.98, 1, 0.99,
       1, 1, 1, 1)

# had to average if both estimates available
X = c(NA, -1.6, NA, 18.3, 
      15.9, NA, NA, 1.0, 
      27.6, 19.65, 24.6, 4.9, 
      NA, 1.68, NA, NA, 
      NA,NA, NA, 19.5,
      NA, -0.7, 18.55, NA,
      NA, 3.6, 5.5, NA,
      16, 1.9, 9, 1.6)

my.data <- data.frame(MT = MT, X = X, UDG = as.factor(UDG))

my.lm <- lm(X ~ MT)
plot(MT, X, xlab = "% Mito. Consensus", ylab = "% X Contam.", main = "R^2 = 0.215, pvalue = 0.052")
abline(my.lm, col = "red", lwd=3)

plot(my.data$UDG, my.data$MT, ylab = "% MT Consensus.")
plot(my.data$UDG, my.data$X, ylab = "% X Contam.")

model <- get(load("../utilities/moderns_Fu/clus_3/model.rda"))
omega <- model$omega[51:94,3]
archaic.names <- rownames(model$omega)[51:94]
matches = lapply(names, function(x) grep(x, archaic.names))
modern.per <- unlist(lapply(matches, function(x) max(omega[x])))
plot(MT, modern.per, col=UDG, pch=20)
abline(lm(modern.per ~ MT))
plot(X, modern.per, col=UDG, pch=20)
abline(lm(modern.per ~ X))
