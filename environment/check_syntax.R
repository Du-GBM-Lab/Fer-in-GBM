script_files <- sort(list.files(
  "scripts",
  pattern = "[.]R$",
  recursive = TRUE,
  full.names = TRUE
))

if (!length(script_files)) {
  stop("No curated R scripts were found under scripts/.")
}

failed <- character()

for (script_file in script_files) {
  result <- tryCatch(
    {
      parse(file = script_file, keep.source = TRUE)
      "PASS"
    },
    error = function(e) {
      failed <<- c(failed, script_file)
      paste("FAIL:", conditionMessage(e))
    }
  )
  cat(sprintf("%-70s %s\n", script_file, result))
}

if (length(failed)) {
  stop("R syntax check failed for: ", paste(failed, collapse = ", "))
}

message("All curated R scripts passed parsing.")
