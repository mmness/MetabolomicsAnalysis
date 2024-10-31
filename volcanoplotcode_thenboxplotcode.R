#Set up data --------------------------------------------------------

##Remember to install install.packages("ggplot2")
library(ggplot2)
library(dplyr)
library(tidyr)
library(openxlsx)
library(officer)
library(rvg)
library(vioplot)


#read in df, enter your number of metadata columns
df<-read.csv("merged_Cohort2_feces.csv", header=TRUE)
nodeTab<-read.csv("CAPCOR-HILICPos-node.csv", header=TRUE)

NAME<-"HILICpos_Cohort2_feces_GSM"
num_metadata_columns<-32

df <- subset(df, Giardiastrain %in% c("GS_M", "No infection")) #TWO groups you want to compare
ToTest<-df$Giardiastrain #edit with your column you want to analyze

obj<-sort( unique(ToTest), decreasing = FALSE ) #change to TRUE if you want to change order

#set significance cutoff
pvalue_sig<-0.05
FC_low<-0.5
FC_high<-2

# TestingNormalization ----------------------------------------------------
f<-num_metadata_columns+1
df_meta <- df[, 1:num_metadata_columns]
df_remaining <- df[, f:ncol(df)]
df_remaining <- df_remaining[, colSums(df_remaining) != 0]
df <- cbind(df_meta, df_remaining)

ncol<- ncol(df)
shapiro_pvals<-numeric(ncol-num_metadata_columns) #number of columns-number of metadata columns 
for (i in 1:length(shapiro_pvals)){
  output<-shapiro.test(df[,i+num_metadata_columns])  #removing the metadata columns.
  shapiro_pvals[i]<-output$p.value}

min(shapiro_pvals)
max(shapiro_pvals)

if (min(shapiro_pvals) < 0.05) {
  print("Data is not normally distributed, use wilcox test")} else {
  print("Data is normally distributed, use t test")}

# FC calc-----------------------------------------------------------------
if(length(obj)>2){print("STOP:column has more than 2 objects :(")}

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
  wilcox_result <- wilcox.test(sub[[col_names[2]]] ~ sub[[col_names[1]]]) #change to t test if normaally distributed
  info_forVP$pval[i] <- wilcox_result$p.value
} 

#adding in some more pval options
info_forVP$pval_adj<- p.adjust(info_forVP$pval, method = "fdr") #can change
info_forVP$negLog10pval<--log10(info_forVP$pval)
info_forVP$negLog10pval_adj<--log10(info_forVP$pval_adj)

# Customize plot ----------------------------------------------------------

###NEED TO EDIT-How to deal with infinities. if you have none or are using the regular FC you can ignore
#Option 1: delete from df
#info_forVP <- info_forVP[info_forVP$FoldChange != 0, ]
#info_forVP <- info_forVP[info_forVP$FoldChange != Inf, ]

#Option 2: transform to a high number and low number. This replaces them with the max and min values but you can edit
max_value <- max(info_forVP$log2FC[info_forVP$log2FC != Inf], na.rm = TRUE)
info_forVP[info_forVP == Inf] <- max_value

min_value <- min(info_forVP$log2FC[info_forVP$log2FC != -Inf], na.rm = TRUE)
info_forVP[info_forVP == -Inf] <- min_value

#choose your axis and cutoffs
y_axis<-info_forVP$negLog10pval #some pvalue variation. Usually -log10pval is used
x_axis<-info_forVP$log2FC #some FC variation. usually log2FC 

SigCuttoff_forFC_low<- log2(0.75)
SigCuttoff_forFC_high<-log2(1.25)
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

ggsave(filename = paste0(NAME, "volcano_plot.png" ), plot = last_plot(), width = 10, height = 7, dpi = 300)

#prints info that clarifies some things
cat("NOTE: on the plot, upregulated means that it is higher in", obj[2], ".", "\n")
cat(sum(data$color == "upregulated"), "are upregulated")
cat(sum(data$color == "downregulated"), "are downregulated")

# Cytoscape-compatible feature information export -----------------------------------
powers<-as.data.frame( cbind(features, means1, means2, sds_1, sds_2, sds, shapiro_pvals) )
colnames(powers)<- c("Feature", paste0(obj[1], "-mean"), paste0(obj[2], "-mean"), paste0(obj[1], "-std dev"), paste0(obj[2], "-std dev"), "std dev all", "shapiro_pval")
feature_info<-merge(info_forVP, powers, by="Feature")

feature_info <- feature_info %>%
  separate(Feature, into = c("mz", "RetentionTime", "ClusterIndex"), sep = "_", remove = FALSE)
feature_info$mz <- as.numeric(sub("X", "", feature_info$mz))
feature_info$RetentionTime <- as.numeric(feature_info$RetentionTime)
feature_info$ClusterIndex <- as.numeric(feature_info$ClusterIndex)
feature_info$'shared name' <- as.numeric(feature_info$ClusterIndex)

write.csv(feature_info, file = paste0(NAME, "_feature_significance_cyto.csv"), row.names = FALSE)

# Boxplots ------------------------------------------------------

#if you want to compare more than two objects- this is not the code for you. 
#using the Feature_info variable created above to ID significant features to analyze

#can change cuttoff here
significant_upreg <- feature_info$Feature[feature_info$pval < pvalue_sig & feature_info$FoldChange > FC_high]
significant_downreg <- feature_info$Feature[feature_info$pval < pvalue_sig & feature_info$FoldChange < FC_low]

allsigs<-c(significant_upreg, significant_downreg)

# plot_list <- list()
# for (i in seq_along(allsigs)) {
#   feat<-allsigs[i]
#   box<-cbind( ToTest,df[feat]  )
#   boxplot(box[[2]] ~ box[[1]],
#           data = box,
#           col = c("lightblue", "lightgreen"),
#           main = feat,
#           xlab = "Group",
#           ylab = "Value",
#           names = unique(box[[1]]),
#           border = "black",
#           notch = FALSE
#           )
#   stripchart(box[[2]] ~ box[[1]],
#              data = box,
#              method = "jitter",
#              jitter = 0.05,
#              pch = 19,
#              col = "darkblue",
#              vertical = TRUE,
#              add = TRUE)
#   plot_list[[i]] <- recordPlot() 
# }

# boxplot-powerpoint export -----------------------------------------------
# ppt <- read_pptx()
# 
# for (i in 1:length(plot_list)) {
#   ppt <- add_slide(ppt, layout = "Title and Content", master = "Office Theme")
#   img_file <- tempfile(fileext = ".png")
#   png(filename = img_file, width = 1000, height = 800, res=300)  #can change resolution
#   replayPlot(plot_list[[i]])  
#   dev.off()
#   ppt <- ph_with(ppt, external_img(img_file), 
#                  location = ph_location(left = 0.5, top = 1, width = 8, height = 6))  
# }
# 
# print(ppt, target = "presentation_allfeatureboxplots.pptx")

# violin plot instead -------------------------------------------------------------

# plot_list <- list()
# for (i in seq_along(allsigs)) {
#   feat<-allsigs[i]
#   box<-cbind( ToTest,df[feat]  )
#   vioplot(box[[2]] ~ box[[1]],
#           col = c("lightcoral", "lightblue3"),  # Colors for each group
#           main = feat,
#           xlab = "Group",
#           ylab = "Value",
#           names = unique(box[[1]]),
#           border = "black",
#           lty = 1,
#           drawRect = TRUE
#           )
#   plot_list[[i]] <- recordPlot() 
# }


# violinPlot-powerpoint export --------------------------------------------


# ppt <- read_pptx()
# 
# for (i in 1:length(plot_list)) {
#   ppt <- add_slide(ppt, layout = "Title and Content", master = "Office Theme")
#   img_file <- tempfile(fileext = ".png")
#   png(filename = img_file, width = 1000, height = 800, res=300)  #can change resolution
#   replayPlot(plot_list[[i]])  
#   dev.off()
#   ppt <- ph_with(ppt, external_img(img_file), 
#                  location = ph_location(left = 0.5, top = 1, width = 8, height = 6))  
# }
# 
# print(ppt, target = "presentation_allviolinplots.pptx")

# Boxplots-only annotations -----------------------------------------------

significant_upreg <- feature_info$Feature[feature_info$pval < pvalue_sig & feature_info$FoldChange > FC_high]
significant_downreg <- feature_info$Feature[feature_info$pval < pvalue_sig & feature_info$FoldChange < FC_low]

feature_info <- feature_info %>%
  mutate(Significant_Cuttoff_parametersfromR = case_when(
    Feature %in% significant_upreg ~ "Upregulated",
    Feature %in% significant_downreg ~ "Downregulated",
    TRUE ~ "NotSignificant"  
  ))

nodeTab_selected <- nodeTab %>%
  select(shared.name, MassDiff, Compound_Name, componentindex)
colnames(nodeTab_selected)<-gsub("shared.name", "ClusterIndex", colnames(nodeTab_selected))

feature_info_new <- merge(feature_info, nodeTab_selected, by = "ClusterIndex")
to_plot <- subset(feature_info_new, feature_info_new$Significant_Cuttoff_parametersfromR != "NotSignificant" & feature_info_new$Compound_Name != "")

plot_list <- list()
for (i in 1:nrow(to_plot)) {
  feat<-to_plot$Feature[i]
  box<-cbind( ToTest,df[feat]  )
  wrapped_title <- paste(strwrap(to_plot$Compound_Name[i], width = 50), collapse = "\n")
  boxplot(box[[2]] ~ box[[1]], 
          data = box, 
          col = c("lightblue", "lightgreen"),  
          main = wrapped_title,
          xlab = "Group", 
          ylab = "Value", 
          names = unique(box[[1]]),  
          border = "black",
          cex.main = 0.7,
          notch = FALSE)  
  stripchart(box[[2]] ~ box[[1]], 
             data = box, 
             method = "jitter", 
             jitter = 0.05,  
             pch = 19,      
             col = "darkblue", 
             vertical = TRUE, 
             add = TRUE) 
  plot_list[[i]] <- recordPlot() 
}


# powerpoint-export-annotations -------------------------------------------

ppt <- read_pptx()

for (i in 1:length(plot_list)) {
  ppt <- add_slide(ppt, layout = "Title and Content", master = "Office Theme")
  img_file <- tempfile(fileext = ".png")
  png(filename = img_file, width = 1000, height = 800, res=200)  #can change resolution
  replayPlot(plot_list[[i]])  
  dev.off()
  ppt <- ph_with(ppt, external_img(img_file), 
                 location = ph_location(left = 0.5, top = 1, width = 8, height = 6))  
}

print(ppt, target = paste0(NAME, "_presentation_annotations.pptx"))



# Boxplots-Only for those with analogs --------------------------------

significant_upreg <- feature_info$Feature[feature_info$pval < pvalue_sig & feature_info$FoldChange > FC_high]
significant_downreg <- feature_info$Feature[feature_info$pval < pvalue_sig & feature_info$FoldChange < FC_low]

feature_info <- feature_info %>%
  mutate(Significant_Cuttoff_parametersfromR = case_when(
    Feature %in% significant_upreg ~ "Upregulated",
    Feature %in% significant_downreg ~ "Downregulated",
    TRUE ~ "NotSignificant"  
  ))

nodeTab_selected <- nodeTab %>%
  select(shared.name, MassDiff, Analog.Compound_Name, componentindex)
colnames(nodeTab_selected)<-gsub("shared.name", "ClusterIndex", colnames(nodeTab_selected))

feature_info_new <- merge(feature_info, nodeTab_selected, by = "ClusterIndex")
to_plot <- subset(feature_info_new, feature_info_new$Significant_Cuttoff_parametersfromR != "NotSignificant" & feature_info_new$Analog.Compound_Name != "")

plot_list <- list()
for (i in 1:nrow(to_plot)) {
  feat<-to_plot$Feature[i]
  box<-cbind( ToTest,df[feat]  )
  wrapped_title <- paste(strwrap(to_plot$Analog.Compound_Name[i], width = 50), collapse = "\n")
  boxplot(box[[2]] ~ box[[1]], 
          data = box, 
          col = c("lightblue", "lightgreen"),  
          main = wrapped_title,
          xlab = "Group", 
          ylab = "Value", 
          names = unique(box[[1]]),  
          border = "black",
          cex.main = 0.7,
          notch = FALSE)  
  stripchart(box[[2]] ~ box[[1]], 
             data = box, 
             method = "jitter", 
             jitter = 0.05,  
             pch = 19,      
             col = "darkblue", 
             vertical = TRUE, 
             add = TRUE) 
  plot_list[[i]] <- recordPlot() 
}


# powerpoint-export-analogs -----------------------------------------------

ppt <- read_pptx()

for (i in 1:length(plot_list)) {
  ppt <- add_slide(ppt, layout = "Title and Content", master = "Office Theme")
  img_file <- tempfile(fileext = ".png")
  png(filename = img_file, width = 1000, height = 800, res=300)  #can change resolution
  replayPlot(plot_list[[i]])  
  dev.off()
  ppt <- ph_with(ppt, external_img(img_file), 
                 location = ph_location(left = 0.5, top = 1, width = 8, height = 6))  
}
print(ppt, target = paste0(NAME, "_presentation_analogs.pptx"))


# ExportData --------------------------------------------------------------
nodeTab_selected <- nodeTab %>%
  select(shared.name, MassDiff, Compound_Name, Analog.Compound_Name, componentindex)
colnames(nodeTab_selected)<-gsub("shared.name", "ClusterIndex", colnames(nodeTab_selected))
feature_info_new <- merge(feature_info, nodeTab_selected, by = "ClusterIndex")
feature_info_onlysigs <- subset(feature_info_new, (FoldChange < FC_low | FoldChange > FC_high) & pval < pvalue_sig)


#write.csv(feature_info_new, "allFeature.csv")
write.csv(feature_info_onlysigs, file = paste0(NAME, "Features_onlysig.csv")) 

