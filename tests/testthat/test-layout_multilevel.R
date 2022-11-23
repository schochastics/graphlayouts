test_that("layout_mulitlevel works", {
  data("multilvl_ex")
  expect_is(layout_as_multilevel(multilvl_ex,type = "all", alpha = 25, beta = 45),"matrix")
  expect_is(layout_as_multilevel(multilvl_ex,type = "separate",
                            FUN1 = layout_as_backbone,
                            FUN2 = layout_with_stress,
                            alpha = 25, beta = 45),"matrix")
  expect_is(layout_as_multilevel(multilvl_ex,
                                 type = "fix2",
                                 FUN2 = layout_with_stress,
                                 alpha = 25, beta = 45
  ),"matrix")
  expect_error(layout_as_multilevel(igraph::make_full_graph(10)))
})
