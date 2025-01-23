test_that("layout_mulitlevel works", {
    testthat::skip_if_not_installed("oaqc")
    data("multilvl_ex")
    expect_is(layout_as_multilevel(multilvl_ex, type = "all", alpha = 25, beta = 45), "matrix")
    expect_is(layout_as_multilevel(multilvl_ex,
        type = "separate",
        FUN1 = layout_as_backbone,
        FUN2 = layout_with_stress,
        alpha = 25, beta = 45
    ), "matrix")
    expect_is(layout_as_multilevel(multilvl_ex,
        type = "fix2",
        FUN2 = layout_with_stress,
        alpha = 25, beta = 45
    ), "matrix")
    expect_is(layout_as_multilevel(multilvl_ex,
        type = "fix1",
        FUN1 = layout_with_stress,
        alpha = 25, beta = 45
    ), "matrix")
    expect_error(layout_as_multilevel(igraph::make_full_graph(10)))

    expect_error(layout_as_multilevel(multilvl_ex, type = "fix1"))
    expect_error(layout_as_multilevel(multilvl_ex, type = "fix2"))
    expect_error(layout_as_multilevel(multilvl_ex, type = "fix3"))

    g <- igraph::delete_vertex_attr(multilvl_ex, "lvl")
    expect_error(layout_as_multilevel(g, type = "fix2"))

    expect_error(layout_as_multilevel(multilvl_ex, type = "separate"))
    expect_error(layout_as_multilevel(multilvl_ex,
        type = "separate",
        FUN1 = layout_as_backbone, params1 = list(wrong = "a"),
        FUN2 = layout_with_stress
    ))
    expect_error(layout_as_multilevel(multilvl_ex,
        type = "separate",
        FUN1 = layout_as_backbone,
        FUN2 = layout_with_stress, params2 = list(wrong = "a")
    ))
    expect_no_error(layout_as_multilevel(multilvl_ex,
        type = "separate",
        FUN1 = layout_as_backbone, params1 = list(keep = 0.3),
        FUN2 = layout_with_stress, params2 = list(bbox = 15)
    ))
})
