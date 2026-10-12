test_that("panel functions run without error on jittered data", {
  popgrid <- load_gz_pop_raster()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(
    fig_study_design_panel_map(contacts_V1_jittered_100m, popgrid, add = TRUE)
  )
  expect_no_error(
    fig_study_design_panel_trajectories(contacts_V1_jittered_100m, popgrid,
                                         seed = 1, add = TRUE)
  )
  expect_no_error(
    fig_study_design_panel_distcdf(contacts_V1_jittered_100m, add = TRUE)
  )
})

test_that("fig_study_design composes all three panels without error", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(
    fig_study_design(contacts_V1_jittered_100m, load_gz_pop_raster(), seed = 1)
  )
})
