test_that("fit_pure_gravity returns a sane single-subset fit", {
  res <- fit_pure_gravity(contacts_V1_jittered_100m, load_gz_pop_raster(),
                           datasubset = "ALL", noRepeats = 1)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_true(is.finite(res$lnlike))
  expect_true(res$point >= 0.1 && res$point <= 6)
  expect_true(res$lower <= res$point && res$point <= res$upper)
})

test_that("table1_pure_gravity returns one row per subset", {
  res <- table1_pure_gravity(contacts_V1_jittered_100m, load_gz_pop_raster(),
                              datasubsets = c("ALL", "ADULTS"), noRepeats = 1)
  expect_equal(nrow(res), 2)
  expect_setequal(res$datasubset, c("ALL", "ADULTS"))
  expect_true(all(is.finite(res$lnlike)))
})
