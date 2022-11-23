test_that("layout_manipulate works", {
  xy <- cbind(c(0,0,0,0),c(1,2,3,4))
  xynew <- layout_mirror(layout_mirror(xy,axis = "horizontal"),axis = "horizontal")
  testthat::expect_true(all(xynew==xy))
  xynew <- layout_mirror(layout_mirror(xy,axis = "vertical"),axis = "vertical")
  testthat::expect_true(all(xynew==xy))

  expect_equal(layout_rotate(xy,180)[,2],c(-1,-2,-3,-4))
})
