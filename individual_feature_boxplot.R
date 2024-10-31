library(ggplot2)
library(ggsignif)


merged<-read.csv("merged_HILICpos_Cohort2_colon.csv")
groups<-merged$Giardiastrain 

sig_cutoff<-0.05

#num_vector<- c("6518", "7406") #Kynur
#num_vector<- c("6925", "24892", "572", "10018", "10378", "11194") #arginine
#num_vector<- c("26547", "26731", "7075", "14234", "10361", "16703") #creatine
#num_vector<- c("12016", "6324", "24660", "19496") #tryptophan
#num_vector<- c("6624", "25058", "17361", "6486") #citruline
#num_vector<- c("29477", "19931", "10738", "2349") #xanthine
#num_vector<- c("8259", "25704") #proline
num_vector<- c("21466", "10293") #Riboflavn

# plot ---------------------------------------------------------------------
set<-0.00005

matching_columns <- colnames(merged)[unlist(sapply(num_vector, function(x) {
  matches <- grep(paste0("_", x, "$"), colnames(merged))
  if (length(matches) > 0) matches else NA
}))]
matching_columns <- matching_columns[!is.na(matching_columns)]

summed <- unlist(matching_columns)
box<-cbind(groups,merged[, summed] )
box<-cbind(box, SumData = rowSums(box[, 2:(length(summed)+1)]))

box$groups <- factor(box$groups, levels = c("No infection", "GS_M", "WB"))
box <- box[order(box$groups), ]

combos<- combn(unique(box$groups), 2)
combos<- data.frame(Group1 = combos[1,], Group2 = combos[2,])

combos$pvalue <- NA
for (i in 1:nrow(combos)){
  df <- box[box$groups %in% unlist(combos[i,1:(ncol(combos)-1)]), ]
  wilcox_test <- wilcox.test(SumData ~ groups, data = df)
  combos$pvalue[i] <- wilcox_test$p.value
  }

sig_combos<- combos[combos$pvalue < sig_cutoff, ]
sig_combos <- sig_combos[ , !(names(sig_combos) %in% "pvalue")]

signif_comparisons <- list()

if (nrow(sig_combos) > 1) {
  for (i in 1:nrow(sig_combos)) {
    signif_comparisons[[i]] <- as.vector(unlist(sig_combos[i, ]))
  }
} 

y_positions <- c( (max(box$SumData)+set), (max(box$SumData)+(2*set)), (max(box$SumData)+(3*set)))  
ggplot(box, aes(x = groups, y = SumData, fill = groups)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(width = 0.2, size = 4, color = "grey2", alpha = 0.5) +  
  theme_minimal() +
  labs(y = "Abundance", title = "Colon") +  
  scale_fill_manual(values = c("green3", "orange", "red")) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank()
  ) +
  geom_signif(comparisons = signif_comparisons, 
              annotations = "*", 
              map_signif_level = FALSE, 
              textsize = 5, 
              y_position = y_positions, 
              tip_length = 0.02)




