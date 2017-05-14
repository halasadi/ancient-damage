

pc_data_frame <- cbind.data.frame(PC1 = c(2,4,5, 10, 12, 12),
                                  PC2 = c(-2, 10, 2, 2, 3, 4),
                                  PC3 = c(23, 10, 12, 1, 2, 4))
labs <- c("a", "a", "a", "b", "b", "b")

graphList <- vector(mode="list");

graphList[[1]] <- ggplot2::qplot(pc_data_frame[,1], pc_data_frame[,2], data=pc_data_frame,
                                 colour=factor(labs))  +
  ggplot2::scale_colour_manual(values = cols, guide = ggplot2::guide_legend(title = "Groups")) +
  ggplot2::xlab(paste0(pcs_to_plot[m])) + ggplot2::ylab(paste0(pcs_to_plot[n]))

graphList[[2]] <- ggplot2::qplot(pc_data_frame[,1], pc_data_frame[,3], data=pc_data_frame,
                                 colour=factor(labs))  +
  ggplot2::scale_colour_manual(values = cols, guide = ggplot2::guide_legend(title = "Groups")) +
  ggplot2::xlab(paste0(pcs_to_plot[m])) + ggplot2::ylab(paste0(pcs_to_plot[n]))

library(ggplot2)
graphList <- list();
l <- 1
m <- 1
for(n in 2:3){
  # a <- ggplot2::qplot(pc_data_frame[,m], pc_data_frame[,n],
  #                     colour=factor(labs))  +
  #   ggplot2::scale_colour_manual(values = cols, guide = ggplot2::guide_legend(title = "Groups")) +
  #   ggplot2::xlab(paste0(pcs_to_plot[1])) + ggplot2::ylab(paste0(pcs_to_plot[n]))

  a <- ggplot(pc_data_frame, aes_string(x = paste0("PC",m), y = paste0("PC", n))) +
    geom_point(aes(colour = factor(labs))) + ggplot2::scale_colour_manual(values = cols, guide = ggplot2::guide_legend(title = "Groups")) +
    ggplot2::xlab(paste0(pcs_to_plot[1])) + ggplot2::ylab(paste0(pcs_to_plot[n]))
  graphList[[l]]  <- a
  l <- l+1
}
graphList[[1]]
graphList[[2]]
