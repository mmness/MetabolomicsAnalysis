#Set up your data --------------------------------------------------------

##install if needed
#install.packages("ggplot2")
#install.packages("dplyr")

library(ggplot2)
library(dplyr)

#read in df, enter your number of metadata columns
df<-read.csv("mergedData.csv", header=TRUE)
num_metadata_columns<-5

df <- subset(df, YOURCOLUMN %in% c("Infected", "No infection")) #TWO groups you want to compare
ToTest<-df$YOURCOLUMN #edit with your column you want to analyze

obj<-sort( unique(ToTest), decreasing = FALSE ) #change to TRUE if you want to change order

#now you can run as-is until line 84

# FC calc-----------------------------------------------------------------
if(length(obj)>2){print("STOP:column has more than 2 objects :(")}
df_meta <- df[, 1:num_metadata_columns]

f<-num_metadata_columns+1
df_remaining <- df[, f:ncol(df)]
df_remaining <- df_remaining[, colSums(df_remaining) != 0]

df <- cbind(df_meta, df_remaining)


features<-colnames(df)[-(1:num_metadata_columns)]
info_forVP <- data.frame(Feature = features, FoldChange = numeric(length(features)))

means1<-numeric()
means2<-numeric()
sds_1<-numeric()
sds_2<-numeric()

for (i in seq_along(features)) {
  feature <- features[i]  #change so that it makes a dataframe. then means gets tested with the df
  
  set<-cbind( ToTest   ,df[feature]  )
  
  mean_value_1 <- mean(set[[feature]][set$ToTest == obj[1]], na.rm = TRUE)
  mean_value_2 <- mean(set[[feature]][set$ToTest == obj[2]], na.rm = TRUE)
  
  means1[i]<-mean(set[[feature]][set$ToTest == obj[1]], na.rm = TRUE)
  means2[i]<-mean(set[[feature]][set$ToTest == obj[2]], na.rm = TRUE)
  sds_1[i]<-sd(set[[feature]][set$ToTest == obj[1]], na.rm = TRUE)
  sds_2[i]<-sd(set[[feature]][set$ToTest == obj[2]], na.rm = TRUE)
  
  sds<-sd(set[[feature]])
  
  FC <- mean_value_2 / mean_value_1
  info_forVP$FoldChange[i] <- FC
}  

#adding log2 and log10FC 
info_forVP$log2FC <- log2(info_forVP$FoldChange)
info_forVP$log10FC <- log10(info_forVP$FoldChange)


# pvalue calculations -----------------------------------------------------

#Calculate p-value (x axis_)
for (i in seq_along(features)) {
  feature <- features[i]
  sub<-cbind(ToTest, df[feature])
  col_names <- names(sub)
  wilcox_result <- wilcox.test(sub[[col_names[2]]] ~ sub[[col_names[1]]])
  info_forVP$pval[i] <- wilcox_result$p.value
} 

#adding in some more pval options
info_forVP$pval_adj<- p.adjust(info_forVP$pval, method = "fdr") #can change
info_forVP$negLog10pval<--log10(info_forVP$pval)
info_forVP$negLog10pval_adj<--log10(info_forVP$pval_adj)

# Customize plot ----------------------------------------------------------

###NEED TO EDIT-How to deal with infinities. if you have none or are using the regular FC you can ignore
#Option 1: delete from df

info_forVP <- info_forVP[info_forVP$FoldChange != 0, ]
info_forVP <- info_forVP[info_forVP$FoldChange != Inf, ]

#Option 2: transform to a high number and low number. This replaces them with the max and min values but you can edit
max_value <- max(info_forVP$log2FC[info_forVP$log2FC != Inf], na.rm = TRUE)
info_forVP[info_forVP == Inf] <- max_value

min_value <- min(info_forVP$log2FC[info_forVP$log2FC != -Inf], na.rm = TRUE)
info_forVP[info_forVP == -Inf] <- min_value

#choose your axis and cutoffs
y_axis<-info_forVP$negLog10pval_adj #some pvalue variation. Usually -log10pval is used
x_axis<-info_forVP$log2FC #some FC variation. usually log2FC 

SigCuttoff_forFC_low<- log2(.5)
SigCuttoff_forFC_high<-log2(1.5)
SigCuttoff_forPval<- -log10(0.05)

# Plot --------------------------------------------------------------------
data <- data.frame(x_axis, y_axis)

data <- data %>%
  mutate(color = case_when(
    y_axis > SigCuttoff_forPval & x_axis > SigCuttoff_forFC_high ~ "upregulated",   #Label for blue
    y_axis > SigCuttoff_forPval & x_axis < SigCuttoff_forFC_low ~ "downregulated",  #Label for red
    TRUE ~ "not_significant"  
  )) %>%
  mutate(color = factor(color, levels = c("upregulated", "downregulated", "not_significant")))  # Ensure correct factor levels

ggplot(data, aes(x = x_axis, y = y_axis, color = color)) +
  geom_point() +  
  geom_hline(yintercept = SigCuttoff_forPval, color = "black", linetype = "dashed") +  
  geom_vline(xintercept = SigCuttoff_forFC_low, color = "black", linetype = "dashed") +  
  geom_vline(xintercept = SigCuttoff_forFC_high, color = "black", linetype = "dashed") +  
  scale_x_continuous(limits = c(min(x_axis, na.rm = TRUE), max(x_axis, na.rm = TRUE))) +  
  scale_y_continuous(limits = c(min(y_axis, na.rm = TRUE), max(y_axis, na.rm = TRUE))) +  
  scale_color_manual(values = c("upregulated" = "blue", "downregulated" = "red", "not significant" = "grey")) +  
  guides(color = guide_legend(title = NULL, override.aes = list(size = 4))) +  
  theme_minimal() +
  labs(x = "Log2FoldChange", y = "neg_log10pval", title = "Volcano Plot") +
  theme(plot.title = element_text(hjust = 0.5))  


#prints info that clarifies some things
cat("NOTE: on the plot, upregulated means that it is higher in", obj[2])
cat(sum(data$color == "upregulated"), "are upregulated")
cat(sum(data$color == "downregulated"), "are downregulated")

# Optional: feature information export -----------------------------------
powers<-as.data.frame( cbind(features, means1, means2, sds_1, sds_2, sds) )
colnames(powers)<- c("Feature", paste0(obj[1], "-mean"), paste0(obj[2], "-mean"), paste0(obj[1], "-std dev"), paste0(obj[2], "-std dev"), "std dev all")
feature_info<-merge(info_forVP, powers, by="Feature")

feature_info <- feature_info %>%
  separate(Feature, into = c("mz", "RetentionTime", "ClusterIndex"), sep = "_", remove = FALSE)
feature_info$mz <- as.numeric(sub("X", "", feature_info$mz))
feature_info$RetentionTime <- as.numeric(feature_info$RetentionTime)
feature_info$ClusterIndex <- as.numeric(feature_info$ClusterIndex)
feature_info$'shared name' <- as.numeric(feature_info$ClusterIndex)

write.csv(feature_info, "feature_significance.csv", row.names = FALSE)














        