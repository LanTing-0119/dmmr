set.seed(123)
mat = matrix(rnorm(80, 2), 8, 10)
mat = rbind(mat, matrix(rnorm(40, -2), 4, 10))
rownames(mat) = paste0("R", 1:12)
colnames(mat) = paste0("C", 1:10)

ha_column1 = HeatmapAnnotation(points = anno_points(rnorm(10)), 
                               annotation_name_side = "left")
ht1 = Heatmap(mat, name = "ht1", km = 2, column_title = "Heatmap 1", 
              top_annotation = ha_column1, row_names_side = "left")

ha_column2 = HeatmapAnnotation(type = c(rep("a", 5), rep("b", 5)),
                               col = list(type = c("a" = "red", "b" = "blue")))
ht2 = Heatmap(mat, name = "ht2", row_title = "Heatmap 2", column_title = "Heatmap 2",
              bottom_annotation = ha_column2, column_km = 2)

ht_list = ht1 + ht2
draw(ht_list, row_title = "Heatmap list", column_title = "Heatmap list")
list_components()
draw(ht1)
