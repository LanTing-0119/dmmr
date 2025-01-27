set.seed(123)
Heatmap(matrix(rnorm(100), 10), name = "mat",
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:3))),
        row_km = 4)
