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

