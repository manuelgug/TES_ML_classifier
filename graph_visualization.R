library(igraph)

dcasted_list <- readRDS("coocurrence_matrices_REAL_DATA.RDS")

#select a pair of infections
concurrency_matrix <- dcasted_list$`29` # 46, 47,27, 29, 19 were classified as R
concurrency_matrix

# Convert to long format (edges)
edge_list <- melt(concurrency_matrix)
colnames(edge_list) <- c("Variant", "TimePoint", "Value")

# Keep only edges with a connection (Value = 1)
edge_list <- edge_list[edge_list$Value == 1, c("Variant", "TimePoint")]

# Create graph from edge list
g <- graph_from_data_frame(edge_list, directed = FALSE)

# Assign node colors: Blue for Variants, Red for Time Points
V(g)$color <- ifelse(V(g)$name %in% edge_list$Variant, "blue", "red")

# Plot the graph with descriptive colors
plot(g, 
     vertex.label = ifelse(V(g)$color == "blue", NA, V(g)$name),  # Hide allele names
     vertex.label.cex = 0.8, 
     vertex.size = 8, 
     vertex.color = V(g)$color,  # Assign colors
     edge.width = 1.5, 
     edge.color = "gray", 
     layout = layout_with_fr)  # Fruchterman-Reingold layout

