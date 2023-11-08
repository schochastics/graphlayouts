test_that("layout_with_pmds works", {
    g <- igraph::make_graph(~ a - -b - -c - -d:e:f)
    expect_is(layout_with_pmds(g, 5), "matrix")
    expect_equal(ncol(layout_with_pmds(g, 5, dim = 3)), 3)
    expect_error(layout_with_pmds(g))
    expect_error(layout_with_pmds(g, 10))
    expect_error(layout_with_pmds(1))
    expect_no_error(layout_with_pmds(g, pivots = 5, weights = rep(4, 5)))

    g <- igraph::graph.full(10) + igraph::graph.full(10)
    expect_error(layout_with_pmds(g, 10))
})

test_that("layout_with_sparse_stress works", {
    g <- igraph::make_graph(~ a - -b - -c - -d:e:f)
    expect_is(layout_with_sparse_stress(g, 5), "matrix")
    expect_error(layout_with_sparse_stress(g))
    expect_error(layout_with_sparse_stress(g, 10))
    expect_error(layout_with_sparse_stress(1))

    g <- igraph::graph.full(10) + igraph::graph.full(10)
    expect_error(layout_with_sparse_stress(g, 10))
    g <- igraph::graph.full(10)
    expect_warning(layout_with_sparse_stress(g, pivots = 5, weights = rep(4, 45)))
})
