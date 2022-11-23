test_that("layout_with_pmds works", {
  g <- igraph::make_graph( ~ a -- b -- c -- d:e:f)
  expect_is(layout_with_pmds(g,5),"matrix")
  expect_error(layout_with_pmds(g))
  expect_error(layout_with_pmds(g,10))
})

test_that("layout_with_sparse_stress works", {
  g <- igraph::make_graph( ~ a -- b -- c -- d:e:f)
  expect_is(layout_with_sparse_stress(g,5),"matrix")
  expect_error(layout_with_sparse_stress(g))
  expect_error(layout_with_sparse_stress(g,10))
})
