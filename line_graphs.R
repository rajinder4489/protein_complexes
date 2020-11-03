library(reshape2)
library(ggplot2)


adj_mat <- list.files(path = "bnstruct/", pattern = "*adjmat*", full.names = T)
expr_files <- list.files(path = "bnstruct/", pattern = "*expr.txt", full.names = T)

am <- gsub("_adjmat.txt", "", adj_mat)
ef <- gsub("_expr.txt", "", expr_files)

ef_sel <- expr_files[ef %in% am]

linegraphs <- list()
linegraphs_smooth <- list()

#for(i in 1:length(ef_sel))
for(i in 1:length(expr_files))
{
#  expr_f <- read.table(file=ef_sel[i], sep="\t", stringsAsFactors = F)
  expr_f <- read.table(file=expr_files[i], sep="\t", stringsAsFactors = F)
  expr_f <- as.data.frame(t(expr_f[-1,]))
  
  colnames(expr_f) <- expr_f[1,]
  expr_f <- expr_f[-1,]
#  cmplx <- gsub("bnstruct/|_.+", "", ef_sel[i])
  cmplx <- gsub("bnstruct/|_.+", "", expr_files[i])
  
  expr_f$tp <- c("0", "2", "8", "24", "72", "168", "240", "336")
  
  melted <- melt(expr_f, id.vars = c("tp"))
  melted$tp <- factor(melted$tp,levels = c("0", "2", "8", "24", "72", "168", "240", "336"))

  melted$value <- as.double(melted$value)
  melted$Type <- "Subunit"
  
  melted[melted$variable == cmplx, "Type"] <- "Complex"
  
  linegraphs[[cmplx]] <- ggplot(melted, aes(x=tp, y=value, group=variable)) + 
    geom_line(aes(color = variable, linetype = Type), size = 1) +
    lims(y = c(0, ceiling(max(melted$value) + 10 * max(melted$value)/100))) + 
    labs(y="FPKM", x = "Time points (in hours)", color = "Subunits/Complex")
  
  ggsave(file=paste("bnstruct/", cmplx, "_lineplot_expr", ".jpeg", sep=""), 
         plot = linegraphs[[cmplx]], device = NULL, path = NULL, scale = 1, 
         width = 10, height = 4, units = "in", dpi = 200, limitsize = TRUE)
  
#Smoothed  
  linegraphs_smooth[[cmplx]] <- ggplot(melted, aes(x=tp, y=value, group=variable)) + 
    geom_point(aes(color = variable)) + 
    lims(y = c(0, ceiling(max(melted$value) + 10 * max(melted$value)/100))) + 
    labs(y="FPKM", x = "Time points (in hours)", color = "Subunits/Complex") + 
    geom_smooth()
#   geom_smooth(aes(color="variable", fill="variable"))
  
  ggsave(file=paste("bnstruct/", cmplx, "_lineplotsmooth_expr", ".jpeg", sep=""), 
         plot = linegraphs_smooth[[cmplx]], device = NULL, path = NULL, scale = 1, 
         width = 10, height = 4, units = "in", dpi = 200, limitsize = TRUE)
  
}


#selected cases

p1 <- linegraphs[["CPX-1780"]] + labs(title = "Laminin-521 complex", subtitle = "CPX-1780", tag = "A")
p2 <- linegraphs[["CPX-5150"]] + labs(title = "AP-2 Adaptor complex, alpha2 variant", subtitle = "CPX-5150", tag = "B")
#p3 <- linegraphs[["CPX-132"]] + labs(title = "BAT3 complex", subtitle = "CPX-132", tag = "C")

p_grid <- grid.arrange(p1, p2, layout_matrix = cbind(c(1, 2)))

ggsave(file=paste("bnstruct/", "merged_lineplot_expr-1780-5150", ".jpeg", sep=""), 
       plot = p_grid, device = NULL, path = NULL, scale = 1, 
       width = 10, height = 6, units = "in", dpi = 300, limitsize = TRUE)
