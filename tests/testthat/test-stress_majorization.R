context("Test layout_with_stress() on connected graphs")

requireNamespace("igraph")

test_that("it works on directed connected graph", {
  g <- igraph::make_graph( ~ a -+ b +-+ c -+ d:e:f)
  expect_silent(
    r <- layout_with_stress(g)
  )
  expect_is(r, "matrix")
})



test_that("it works on undirected connected graph", {
  g <- igraph::make_graph( ~ a -- b -- c -- d:e:f)
  expect_silent(
    r <- layout_with_stress(g)
  )
  expect_is(r, "matrix")
})

test_that("it works on undirected connected weighted graph", {
  g <- igraph::make_graph( ~ a -- b -- c -- d:e:f)
  igraph::E(g)$weight <- c(1,2,3,4,5)
  expect_silent(
    r <- layout_with_stress(g,weights = NULL)
  )
  expect_is(r, "matrix")
})


context("layout_with_stress() works on disconnected graphs")

test_that("it works on an isolate", {
  g <- igraph::make_graph( ~ a)
  expect_silent(
    r <- layout_with_stress(g)
  )
  expect_is(r, "matrix")
})

test_that("it works on a graph of 5 isolates", {
  g <- igraph::make_graph( ~ a, b, c, d, e)
  expect_silent(
    r <- layout_with_stress(g)
  )
  expect_is(r, "matrix")
})

test_that("it works on an undirected graph of two connected dyads", {
  g <- igraph::make_graph( ~ a -- b, c -- d)
  expect_silent(
    r <- layout_with_stress(g)
  )
  expect_is(r, "matrix")
})


test_that("it works on an undirected graph of two connected dyads with 5 isolates", {
  g <- igraph::make_graph( ~ a -- b, c -- d, e, f, g, h, i)
  expect_silent(
    r <- layout_with_stress(g)
  )
  expect_is(r, "matrix")
})


context("Test layout_with_stress3D() on connected graphs")

test_that("it works on directed connected graph", {
  g <- igraph::make_graph( ~ a -+ b +-+ c -+ d:e:f)
  expect_silent(
    r <- layout_with_stress3D(g)
  )
  expect_is(r, "matrix")
})

test_that("it works on undirected connected weighted graph", {
  g <- igraph::make_graph( ~ a -- b -- c -- d:e:f)
  igraph::E(g)$weight <- c(1,2,3,4,5)
  expect_silent(
    r <- layout_with_stress3D(g,weights = NULL)
  )
  expect_is(r, "matrix")
})

test_that("it works on undirected connected graph", {
  g <- igraph::make_graph( ~ a -- b -- c -- d:e:f)
  expect_silent(
    r <- layout_with_stress3D(g)
  )
  expect_is(r, "matrix")
})


context("layout_with_stress3D() works on disconnected graphs")

test_that("it works on an isolate", {
  g <- igraph::make_graph( ~ a)
  expect_silent(
    r <- layout_with_stress(g)
  )
  expect_is(r, "matrix")
})

test_that("it works on a graph of 5 isolates", {
  g <- igraph::make_graph( ~ a, b, c, d, e)
  expect_silent(
    r <- layout_with_stress3D(g)
  )
  expect_is(r, "matrix")
})

test_that("it works on an undirected graph of two connected dyads", {
  g <- igraph::make_graph( ~ a -- b, c -- d)
  expect_silent(
    r <- layout_with_stress3D(g)
  )
  expect_is(r, "matrix")
})


test_that("it works on an undirected graph of two connected dyads with 5 isolates", {
  g <- igraph::make_graph( ~ a -- b, c -- d, e, f, g, h, i)
  expect_silent(
    r <- layout_with_stress3D(g)
  )
  expect_is(r, "matrix")
})

context("Test layout_with_focus() on connected graphs")

test_that("it works on undirected connected graphs",{
  g <- igraph::make_star(5,mode = "undirected",center = 1)
  expect_error(layout_with_focus(g))
  expect_silent(
    r <- layout_with_focus(g,v = 1)$xy
  )
  expect_is(r, "matrix")
})


test_that("it fails for disconnected graphs",{
  expect_error(layout_with_focus(igraph::graph.empty(n=10,directed = FALSE)))
})

context("Test layout_with_centrality() on connected graphs")

test_that("it works on undirected connected graphs",{
  g <- igraph::make_star(5,mode = "undirected",center = 1)
  expect_error(layout_with_centrality(g))
  expect_silent(
    r <- layout_with_centrality(g,cent = igraph::degree(g))
  )
  expect_is(r, "matrix")
})


test_that("it fails for disconnected graphs",{
  expect_error(layout_with_centrality(igraph::graph.empty(n=10,directed = FALSE)))
})

context("Test layout_with_constrained_stress() on connected graphs")


test_that("it works on undirected connected graph", {
  g <- igraph::make_graph( ~ a -- b -- c -- d:e:f)
  expect_error(layout_with_constrained_stress(g))
  expect_silent(
    r <- layout_with_constrained_stress(g,coord=rep(1,6))
  )
  expect_is(r, "matrix")
})


context("Test layout_with_constrained_stress3D() on connected graphs")
test_that("it works on undirected connected graph", {
  g <- igraph::make_graph( ~ a -- b -- c -- d:e:f)
  expect_error(layout_with_constrained_stress3D(g))
  expect_silent(
    r <- layout_with_constrained_stress3D(g,coord=rep(1,6))
  )
  expect_is(r, "matrix")
})
