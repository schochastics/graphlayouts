test_that("adjust_dist works", {
  D <- matrix(c(0,1,1,Inf,Inf,
                1,0,1,Inf,Inf,
                1,1,0,Inf,Inf,
                Inf,Inf,Inf,0,1,
                Inf,Inf,Inf,1,0),5,5,byrow=TRUE)
  Dlist <- list(D,D)
  expect_no_error(adjust_dist(Dlist))
})

