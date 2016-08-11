

#######   Factor Analysis on Lindo2016 data  ####################

signature_counts <- get(load("../summary_data/signature-counts-clubbed-Lindo2016.rda"))

voom_signature_counts <- t(limma::voom(t(signature_counts))$E);
voom_weights <- t(limma::voom(t(signature_counts))$weights);

Rcpp::sourceCpp("../../benchmarking/flash/cloglik.cpp")
Rcpp::sourceCpp(("../../benchmarking/flash/ATM.cpp"))
source("../../benchmarking/flash/flash.R")
source("../../benchmarking/flash/greedy.R")

ll <- greedy(voom_signature_counts, K=10, 
             flash_para = list(tol=1e-3, maxiter_r1 = 100,
             partype="known", sigmae2_true = voom_weights,
             nonnegative=FALSE));

ll2 <- flashpool(voom_signature_counts, K=6, 
                 tol=1e-3, maxiter_r1 = 100,
                 partype="known", sigmae2_true = voom_weights,
                 nonnegative=FALSE)

### greedy flash finds 6 factors 

# save(ll, file="../rda/greedy_flash_voom_signature_Lindo2016.rda")
# save(ll2, file="../rda/flashpool_flash_voom_signature_Lindo2016.rda")

ll <- get(load("../rda/greedy_flash_voom_signature_Lindo2016.rda"))
ll2 <- get(load("../rda/flashpool_flash_voom_signature_Lindo2016.rda"))

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(ll$l))),
  label = factor(c(rep("Ancient",25), rep("Modern",25)))
)

rownames(ll$l) <- paste0("X", c(1:NROW(ll$l)));

FactorGGBar(loadings = ll$l,
             annotation = annotation,
            palette = list("mid"="white", 
                           "low"="red", 
                           "high"="blue", 
                           "midpoint"=0),
             yaxis_label = "Population Type",
             figure_title = " ",
             axis_tick = list(axis_ticks_length = .1,
                              axis_ticks_lwd_y = .1,
                              axis_ticks_lwd_x = .1,
                              axis_label_size = 7,
                              axis_label_face = "bold"),
            panel=list(panel_rows=2,
                       panel_title="greedy FLASH Factor Loadings Bar plot (with C to T)",
                       panel_title_fontsize=10,
                       panel_title_font=3))

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(ll2$l))),
  label = factor(c(rep("Ancient",25), rep("Modern",25)))
)

rownames(ll2$l) <- paste0("X", c(1:NROW(ll2$l)));

FactorGGBar(loadings = ll2$l,
            annotation = annotation,
            palette = list("mid"="white", 
                           "low"="red", 
                           "high"="blue", 
                           "midpoint"=0),
            yaxis_label = "Population Type",
            figure_title = " ",
            axis_tick = list(axis_ticks_length = .1,
                             axis_ticks_lwd_y = .1,
                             axis_ticks_lwd_x = .1,
                             axis_label_size = 7,
                             axis_label_face = "bold"),
            panel=list(panel_rows=2,
                       panel_title="pool FLASH Factor Loadings Bar plot (with C to T)",
                       panel_title_fontsize=10,
                       panel_title_font=3))

#######  Factor analysis after removing C -> T ##################

signature_counts <- get(load("../summary_data/signature-counts-clubbed-Lindo2016.rda"))

sig_list <- colnames(signature_counts);
sig_list_split <- do.call(rbind, lapply(colnames(signature_counts), function(x) strsplit(as.character(x), split="")[[1]]))

indices <- which(sig_list_split[,3]=="C" & sig_list_split[,6]=="T")

signature_counts_no_CtoT <- signature_counts[,-indices];

voom_signature_counts_no_CtoT <- t(limma::voom(t(signature_counts_no_CtoT))$E);
voom_weights <- t(limma::voom(t(signature_counts_no_CtoT))$weights);

ll <- greedy(voom_signature_counts_no_CtoT, K=10, 
             flash_para = list(tol=1e-3, maxiter_r1 = 100,
                               partype="known", sigmae2_true = voom_weights,
                               nonnegative=FALSE));

ll2 <- flashpool(voom_signature_counts_no_CtoT, K=6, 
                 tol=1e-3, maxiter_r1 = 100,
                 partype="known", sigmae2_true = voom_weights,
                 nonnegative=FALSE)

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(ll$l))),
  label = factor(c(rep("Ancient",25), rep("Modern",25)))
)

rownames(ll$l) <- paste0("X", c(1:NROW(ll$l)));

FactorGGBar(loadings = ll$l,
            annotation = annotation,
            palette = list("mid"="white", 
                           "low"="red", 
                           "high"="blue", 
                           "midpoint"=0),
            yaxis_label = "Population Type",
            figure_title = " ",
            axis_tick = list(axis_ticks_length = .1,
                             axis_ticks_lwd_y = .1,
                             axis_ticks_lwd_x = .1,
                             axis_label_size = 7,
                             axis_label_face = "bold"),
            panel=list(panel_rows=2,
                       panel_title="greedy FLASH Factor Loadings Bar plot (without C to T)",
                       panel_title_fontsize=10,
                       panel_title_font=3))

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(ll2$l))),
  label = factor(c(rep("Ancient",25), rep("Modern",25)))
)

rownames(ll2$l) <- paste0("X", c(1:NROW(ll2$l)));

FactorGGBar(loadings = ll2$l,
            annotation = annotation,
            palette = list("mid"="white", 
                           "low"="red", 
                           "high"="blue", 
                           "midpoint"=0),
            yaxis_label = "Population Type",
            figure_title = " ",
            axis_tick = list(axis_ticks_length = .1,
                             axis_ticks_lwd_y = .1,
                             axis_ticks_lwd_x = .1,
                             axis_label_size = 7,
                             axis_label_face = "bold"),
            panel=list(panel_rows=2,
                       panel_title="pool FLASH Factor Loadings Bar plot (without C to T)",
                       panel_title_fontsize=10,
                       panel_title_font=3))

################  SFA on Lindo 2016 data    ##################

write.table(voom_signature_counts, file="../summary_data/sfa_input_voom_sig_with_CtoT_lindo2016.txt", 
            row.names=FALSE,col.names=FALSE, quote=FALSE, sep="\t")

write.table(voom_signature_counts_no_CtoT, file="../summary_data/sfa_input_voom_sig_without_CtoT_lindo2016.txt", 
            row.names=FALSE,col.names=FALSE, quote=FALSE, sep="\t")

# ./sfa_mac -gen ../../sfa_inputs/sfa_input_voom_sig_without_CtoT_lindo2016.txt -g 50 -n 1536 -k 5 -iter 500 -r 800 -mn -mg -o ../../sfa_outputs/Lindo2016/without_CtoT/voom_lindo_without_CtoT
#  ./sfa_mac -gen ../../sfa_inputs/sfa_input_voom_sig_with_CtoT_lindo2016.txt -g 50 -n 1536 -k 5 -iter 500 -r 800 -mn -mg -o ../../sfa_outputs/Lindo2016/voom_lindo_with_CtoT

###  Case 1 : with CtoT included #######################
# Expected complete Log likelihood at iteration 500: 8047.85
# Marginal log likelihood at iteration 500: -25484
# Residual variance at iteration 500: 72.649
# Residual sum of squares at iteration 500: 7227.68

###   Case 2 : without C to T  ########################
# Expected complete Log likelihood at iteration 500: 12159.9
# Marginal log likelihood at iteration 500: -22889
# Residual variance at iteration 500: 71.2059
# Residual sum of squares at iteration 500: 6649.5

ll_load <- read.table("../utilities/sfa_Lindo2016/with_CtoT/voom_lindo_with_CtoT_lambda.out")
dim(ll_load)

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(ll_load))),
  label = factor(c(rep("Ancient",25), rep("Modern",25)))
)

rownames(ll_load) <- paste0("X", c(1:NROW(ll$l)));

FactorGGBar(loadings = ll_load,
            annotation = annotation,
            palette = list("mid"="white", 
                           "low"="red", 
                           "high"="blue", 
                           "midpoint"=0),
            yaxis_label = "Population Type",
            figure_title = " ",
            axis_tick = list(axis_ticks_length = .1,
                             axis_ticks_lwd_y = .1,
                             axis_ticks_lwd_x = .1,
                             axis_label_size = 7,
                             axis_label_face = "bold"),
            panel=list(panel_rows=2,
                       panel_title="SFA Factor Loadings Bar plot (with C to T)",
                       panel_title_fontsize=10,
                       panel_title_font=3))


ll2_load <- read.table("../utilities/sfa_Lindo2016/without_CtoT/voom_lindo_without_CtoT_lambda.out")
dim(ll2_load)

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(ll2_load))),
  label = factor(c(rep("Ancient",25), rep("Modern",25)))
)

rownames(ll2_load) <- paste0("X", c(1:NROW(ll2_load)));

FactorGGBar(loadings = ll2_load,
            annotation = annotation,
            palette = list("mid"="white", 
                           "low"="red", 
                           "high"="blue", 
                           "midpoint"=0),
            yaxis_label = "Population Type",
            figure_title = " ",
            axis_tick = list(axis_ticks_length = .1,
                             axis_ticks_lwd_y = .1,
                             axis_ticks_lwd_x = .1,
                             axis_label_size = 7,
                             axis_label_face = "bold"),
            panel=list(panel_rows=2,
                       panel_title="SFA Factor Loadings Bar plot (without C to T)",
                       panel_title_fontsize=10,
                       panel_title_font=3))