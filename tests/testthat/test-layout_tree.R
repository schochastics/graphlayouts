test_that("equalangle works on a balanced tree", {
  g <- igraph::make_tree(15, 2, mode = "undirected")
  expect_silent(
    r <- layout_as_tree_unrooted(g)
  )
  expect_is(r, "matrix")
  expect_equal(dim(r), c(15, 2))
  expect_true(all(is.finite(r)))
})

test_that("all modes work and return n x 2 finite matrices", {
  g <- igraph::make_tree(20, 3, mode = "undirected")
  for (m in c("equalangle", "equaldaylight", "stress")) {
    expect_silent(
      r <- layout_as_tree_unrooted(g, mode = m)
    )
    expect_is(r, "matrix")
    expect_equal(dim(r), c(20, 2))
    expect_true(all(is.finite(r)))
  }
})

test_that("it works on a star", {
  g <- igraph::make_star(8, mode = "undirected", center = 1)
  expect_silent(
    r <- layout_as_tree_unrooted(g, mode = "equaldaylight")
  )
  expect_is(r, "matrix")
  expect_equal(dim(r), c(8, 2))
})

test_that("it works on an unbalanced caterpillar tree", {
  g <- igraph::make_graph(~ 1 - 2, 2 - 3, 3 - 4, 4 - 5, 2 - 6, 3 - 7, 4 - 8)
  expect_silent(
    r <- layout_as_tree_unrooted(g, mode = "equaldaylight")
  )
  expect_is(r, "matrix")
  expect_equal(dim(r), c(igraph::vcount(g), 2))
  expect_true(all(is.finite(r)))
})

test_that("branch lengths are respected", {
  g <- igraph::make_star(5, mode = "undirected", center = 1)
  igraph::E(g)$weight <- c(1, 2, 3, 4)
  r <- layout_as_tree_unrooted(g, weights = NULL)
  # distance from center (node 1) to each leaf should match the branch length
  d <- sqrt(rowSums((r[-1, ] - matrix(r[1, ], 4, 2, byrow = TRUE))^2))
  expect_equal(sort(d), c(1, 2, 3, 4))
})

test_that("weights = NA ignores the weight attribute", {
  g <- igraph::make_tree(10, 2, mode = "undirected")
  igraph::E(g)$weight <- runif(igraph::ecount(g), 1, 5)
  g2 <- igraph::delete_edge_attr(g, "weight")
  expect_equal(
    layout_as_tree_unrooted(g, weights = NA),
    layout_as_tree_unrooted(g2, weights = NA)
  )
  expect_false(isTRUE(all.equal(
    layout_as_tree_unrooted(g, weights = NA),
    layout_as_tree_unrooted(g, weights = NULL)
  )))
})

test_that("it works on a forest", {
  g <- igraph::disjoint_union(
    igraph::make_tree(7, 2, mode = "undirected"),
    igraph::make_star(5, mode = "undirected")
  )
  expect_silent(
    r <- layout_as_tree_unrooted(g)
  )
  expect_is(r, "matrix")
  expect_equal(dim(r), c(12, 2))
  expect_true(all(is.finite(r)))
})

test_that("small trees work", {
  expect_equal(dim(layout_as_tree_unrooted(igraph::make_graph(~a))), c(1, 2))
  expect_equal(dim(layout_as_tree_unrooted(igraph::make_graph(~ a - -b))), c(2, 2))
})

test_that("errors on invalid input", {
  expect_error(layout_as_tree_unrooted(5))
  # a graph with a cycle is not a tree
  expect_error(layout_as_tree_unrooted(igraph::make_ring(5)))
  expect_error(layout_as_tree_unrooted(igraph::make_full_graph(4)))
})

test_that("igraph wrapper returns a data frame with coordinates", {
  g <- igraph::make_tree(10, 2, mode = "undirected")
  igraph::V(g)$name <- letters[1:10]
  nodes <- layout_igraph_tree_unrooted(g)
  expect_s3_class(nodes, "data.frame")
  expect_true(all(c("x", "y", "circular") %in% names(nodes)))
  expect_equal(nrow(nodes), 10)
})
