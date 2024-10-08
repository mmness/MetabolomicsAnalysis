FeatureTab<-read.table("QiimeFormat_feattable.txt", sep="\t", header = TRUE)
MetaDat<-read.table("QiimeFormat_metadata.txt", sep = "\t", header = TRUE)

library(vegan)
library(plotly)

sample_names <- FeatureTab[, 1]  
data_matrix <- FeatureTab[, -1]  
MetaDat$SampleID <- MetaDat[, 1]

pcoa_meta <- merge(data.frame(SampleID = sample_names), MetaDat, by = "SampleID")
dist_matrix <- vegdist(data_matrix, method = "bray")
pcoa_results <- cmdscale(dist_matrix, eig = TRUE, k = 3)

#Create a data frame. change to your metadatacolumn
pcoa_df <- data.frame(SampleID = pcoa_meta$SampleID,
                      PC1 = pcoa_results$points[, 1],
                      PC2 = pcoa_results$points[, 2],
                      PC3 = pcoa_results$points[, 3],
                      Tissue = pcoa_meta$Type)

#Calculate variance 
var_explained <- round(100 * pcoa_results$eig / sum(pcoa_results$eig), 1)

#Plot. change parameters if needed
plot_ly(pcoa_df, 
        x = ~PC1, y = ~PC2, z = ~PC3, 
        color = ~Tissue,  
        text = ~SampleID,  
        type = "scatter3d", mode = "markers") %>%
  layout(title = "3D PCoA Plot by Organ",
         scene = list(xaxis = list(title = paste0("PC1 (", var_explained[1], "%)")),
                      yaxis = list(title = paste0("PC2 (", var_explained[2], "%)")),
                      zaxis = list(title = paste0("PC3 (", var_explained[3], "%)"))))
















