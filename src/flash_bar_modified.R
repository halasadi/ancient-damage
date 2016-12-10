

###  FLASH graphs with PVE

signature_counts <- get(load("../summary_data/signature-counts-clubbed-Lindo2016.rda"))

sig_list <- colnames(signature_counts);
sig_list_split <- do.call(rbind, lapply(colnames(signature_counts), function(x) strsplit(as.character(x), split="")[[1]]))

indices <- which(sig_list_split[,3]=="C" & sig_list_split[,6]=="T")

signature_counts_no_CtoT <- signature_counts[,-indices];

voom_signature_counts_no_CtoT <- t(limma::voom(t(signature_counts_no_CtoT))$E);
voom_weights <- t(limma::voom(t(signature_counts_no_CtoT))$weights);

### The FLASH with partype="constant" works better

ll <- greedy(voom_signature_counts_no_CtoT, K=10, 
             flash_para = list(tol=1e-3, maxiter_r1 = 100,
                               partype="constant",
                               nonnegative=FALSE));

ll2 <- flashpool(voom_signature_counts_no_CtoT, K=10, 
                 tol=1e-3, maxiter_r1 = 100,
                 partype="constant",
                 nonnegative=FALSE)


flash_factor_postprocess <- function(ll, data){
  PVE <- array(0, (dim(ll$f)[2]-1));
  sparsity <- array(0, dim(ll$f)[2])
  big_genes_2 <- array(0, dim(ll$f)[2])
  big_genes_5 <- array(0, dim(ll$f)[2])
  prop_positive <- array(0, dim(ll$f)[2])
  prop_negative <- array(0, dim(ll$f)[2])
  for(num in 1:dim(ll$f)[2]){
    max_gene <- max(ll$f[,num]);
    abs_fac <- abs(ll$f[,num]);
    abs_fac_max <- max(abs_fac);
    zero_indices <- which(((abs_fac)/abs_fac_max) < 1e-04);
    ll$f[zero_indices,num]=0;
    abs_fac[zero_indices]=0;
    big_genes_2[num] <- length(which(abs_fac > (0.2*abs_fac_max)))/(length(abs_fac));
    big_genes_5[num] <- length(which(abs_fac > (0.5*abs_fac_max)))/(length(abs_fac));
    prop_positive[num] <- length(which(ll$f[,num] > 0))/length(ll$f[,num]);
    prop_negative[num] <- length(which(ll$f[,num] < 0))/length(ll$f[,num]);
    sparsity[num] <- (length(which(abs_fac==0)))/(length(abs_fac))
    if(num==1){
      temp <- ll$l[,num]%*%t(ll$f[,num]);
      res <- data - temp;
    }else{
      temp <- ll$l[,num]%*%t(ll$f[,num]);
      PVE[num-1] <- (sum(temp^2))/(sum(res^2));
    }
  }
  
  post <- list("PVE"=PVE,
               "sparsity_prop"=sparsity,
               "prop_positive_features"=prop_positive,
               "prop_negative_features"=prop_negative,
               "big_genes_2"=big_genes_2,
               "big_genes_5"=big_genes_5)
  return(post)
}


postprocess_ll <- flash_factor_postprocess(ll, voom_signature_counts_no_CtoT)
pve_percentage <- c("", paste0(": PVE : ", floor(postprocess_ll$PVE*100), "%"))

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(ll$l))),
  label = factor(c(rep("Ancient",25), rep("Modern",25)))
)

rownames(ll$l) <- paste0("X", c(1:NROW(ll$l)));

library(ggplot2)
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
            legend_labels = NULL,
            scale=TRUE,
            panel=list(panel_rows=2,
                       panel_title="greedy FLASH Factor Loadings Bar plot (without C to T)",
                       panel_title_fontsize=10,
                       panel_title_font=3))

scale_ll <- apply(ll$l,2,function(x) x/sd(x));

FactorGGplot(loadings = ll$l[,-1],
             annotation = annotation,
             palette = c( RColorBrewer::brewer.pal(8, "Accent"),
                          RColorBrewer::brewer.pal(4, "Spectral")),
             figure_title = "greedy FLASH Factor loadings stacked barchart",
             yaxis_label = "Factor type",
             order_sample = TRUE,
             sample_order_decreasing = TRUE,
             split_line = list(split_lwd = 0.2,
                               split_col = "black"),
             plot_labels = TRUE,
             legend_labels = pve_percentage[-1],
             scale=TRUE,
             axis_tick = list(axis_ticks_length = .1,
                              axis_ticks_lwd_y = .1,
                              axis_ticks_lwd_x = .1,
                              axis_label_size = 3,
                              axis_label_face = "bold"))
