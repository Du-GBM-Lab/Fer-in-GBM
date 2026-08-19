# Copy this file to config/paths.local.R and edit only the local copy.
# The local file is ignored by Git.

project_paths <- list(
  data_root = "D:/replace/with/local/data",
  results_root = "D:/replace/with/local/results",
  gbmap_reference = "D:/replace/with/GBmap/reference.rds",
  human_spatial_root = "D:/replace/with/human/spatial/data",
  gliomagenesis_object = "D:/replace/with/gliomagenesis/object.rds",
  snatac_object = "D:/replace/with/Combined_snATAC_Annotated_AddMotifsed_RunChromVARed.rds"
)

required_path_names <- names(project_paths)
missing_placeholders <- vapply(
  project_paths,
  function(x) grepl("replace/with", x, fixed = TRUE),
  logical(1)
)

if (any(missing_placeholders)) {
  warning(
    "Replace placeholder paths for: ",
    paste(required_path_names[missing_placeholders], collapse = ", ")
  )
}
