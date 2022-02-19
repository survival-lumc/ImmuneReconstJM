# Targets after DLI - update: in ITT subset

# To-do:
# - correlated long sub-models
# - no more need for predicting from previous model
# - use same specification as pre-DLI (see new proposal)

postDLI_targets <- list(
  #tar_target(postDLI_CD3_long, run_preDLI_longitudinal(...)),
  # For CD4 and CD8 targets..

  # Post-DLI cox: single model needed
  #tar_target(postDLI_cox, run_postDLI_cox()),

  # Joint models with only current val (try CD3 with slope)
  # .. later with jmbayes2 tests

  # Also run tests with gvhd a week earlier..
)
