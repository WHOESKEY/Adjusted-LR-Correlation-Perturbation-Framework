# Reassemble the archived public TPM matrices and retrieve the versioned LR reference.
# Run from the repository root. Generated inputs overwrite results/00_Input only.
library(qs2)
dir.create("results/00_Input", recursive = TRUE, showWarnings = FALSE)
sample_map <- qs_read("data/tcga_selected_sample_map.qs2")
for (cohort in c("GTEx", "TCGA")) {
  files <- list.files(file.path("data/tpm", cohort), pattern = "^samples_[0-9]+\\.qs2$",
                      full.names = TRUE)
  blocks <- lapply(files, qs_read)
  stopifnot(length(blocks) > 0, all(vapply(blocks, function(x) {
    identical(rownames(x), rownames(blocks[[1]]))
  }, logical(1))))
  tpm <- do.call(cbind, blocks)
  stopifnot(!anyDuplicated(colnames(tpm)), !anyDuplicated(rownames(tpm)))
  if (cohort == "GTEx") {
    stopifnot(identical(dim(tpm), c(57853L, 898L)))
    qs_save(tpm, "results/00_Input/gtex_colon_expression_tpm.qs2")
  } else {
    stopifnot(identical(dim(tpm), c(38624L, 668L)), setequal(colnames(tpm), sample_map$sample_id))
    qs_save(tpm, "results/00_Input/tcga_selected_expression_tpm.qs2")
  }
}
# The upstream snapshot includes multiple resources. Step 02 uses 17 human collections.
target <- "results/00_Input/omni_resources.rds"
digest <- "b807279afdb1b2f34fa83a3dd51bbcb5"
installed <- system.file("omni_resources.rds", package = "liana")
if (!file.exists(target) || unname(tools::md5sum(target)) != digest) {
  if (nzchar(installed) && unname(tools::md5sum(installed)) == digest) {
    stopifnot(file.copy(installed, target, overwrite = TRUE))
  } else {
    url <- paste0("https://raw.githubusercontent.com/saezlab/liana/",
                  "6cab46c54234f861ea176c3de77c4b8aa45ecb3d/inst/omni_resources.rds")
    download.file(url, target, mode = "wb")
  }
}
stopifnot(identical(unname(tools::md5sum(target)), digest))
