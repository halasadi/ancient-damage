

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
                       panel_title="Factor Loadings Bar plot",
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
                       panel_title="Factor Loadings Bar plot",
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
                       panel_title="Factor Loadings Bar plot",
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
                       panel_title="Factor Loadings Bar plot",
                       panel_title_fontsize=10,
                       panel_title_font=3))
