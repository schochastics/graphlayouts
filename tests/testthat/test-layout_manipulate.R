test_that("layout_manipulate works", {
  xy <- cbind(c(0,0,0,0),c(1,2,3,4))
  xynew <- layout_mirror(layout_mirror(xy,axis = "horizontal"),axis = "horizontal")
  expect_true(all(xynew==xy))
  xynew <- layout_mirror(layout_mirror(xy,axis = "vertical"),axis = "vertical")
  expect_true(all(xynew==xy))

  expect_error(layout_mirror(xy,"diagonal"))
  expect_error(layout_mirror(c(1,1),"horizontal"))

  expect_equal(layout_rotate(xy,180)[,2],c(-1,-2,-3,-4))
  expect_error(layout_rotate(c(1,1),70))
  expect_error(layout_rotate(xy,"70"))
})
