test_that("reorder_edges works", {
   g <- igraph::make_full_graph(3)
   igraph::E(g)$weight <- c(1,2,3)
   g1 <- graphlayouts::reorder_edges(g,"weight")
   expect_equal(igraph::E(g)$weight,rev(igraph::E(g1)$weight))
})
