#' Release metadata of the fitted model(s) shipped with this package
#'
#' Reads `inst/model/weights.json`: one entry per model variant with the release
#' tag, the download URL of the parameter file, its sha256 and the FINN commit it
#' was fitted under. The first entry is the default variant.
#'
#' @return A named list, one element per variant.
#' @export
model_info = function() {
  path = system.file("model", "weights.json", package = utils::packageName())
  jsonlite::read_json(path, simplifyVector = FALSE)$variants
}

#' Species of the model
#'
#' @param variant Model variant (see [model_info()]); `NULL` takes the default (first).
#' @return A data.frame with the species index used by that variant and its name.
#' @export
species_table = function(variant = NULL) {
  if (is.null(variant)) variant = names(model_info())[1]
  sp = utils::read.csv(system.file("model", "species.csv", package = utils::packageName()))
  sp = sp[sp$variant == variant, c("species", "name")]
  rownames(sp) = NULL
  sp
}

#' Load the fitted model
#'
#' Downloads the fitted parameters of the release this package version belongs
#' to (once; cached in the user cache directory), checks their sha256, and
#' rebuilds the model with the installed FINN through [build_model()]. The model
#' is rebuilt rather than restored with `torch_load()` of a whole model object,
#' because a restored object keeps the methods of the FINN version it was saved
#' under.
#'
#' @param variant Name of the model variant; `NULL` takes the default (first).
#' @param path Optional local parameter file; skips the download.
#' @return A fitted `FINN` model ready for [FINN::simulateForest()].
#' @export
load_model = function(variant = NULL, path = NULL) {
  info = model_info()
  if (is.null(variant)) variant = names(info)[1]
  if (!variant %in% names(info))
    stop("Unknown variant '", variant, "'; available: ", paste(names(info), collapse = ", "), call. = FALSE)
  v = info[[variant]]
  if (is.null(path)) {
    path = file.path(tools::R_user_dir(utils::packageName(), "cache"), v$tag, v$file)
    if (!file.exists(path)) {
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      utils::download.file(v$url, path, mode = "wb")
    }
  }
  sha = digest::digest(file = path, algo = "sha256")
  if (!identical(sha, v$sha256))
    stop("sha256 of ", path, " does not match release ", v$tag, " (", sha, " vs ", v$sha256, ").", call. = FALSE)
  build_model(variant = variant, parameters = torch::torch_load(path))
}

#' Simulate forest dynamics with the fitted model
#'
#' A thin wrapper around [FINN::simulateForest()] with the settings the model
#' was fitted with (patch size, number of patches) as defaults.
#'
#' @param env Environment table with the covariates listed in the model card
#'   (see `vignette("quickstart")`).
#' @param init_cohort Initial cohorts from [FINN::makeInitCohorts()]; `NULL`
#'   starts from bare ground.
#' @param variant Model variant, see [model_info()].
#' @param batchsize Sites per simulation batch; the default is the value pinned
#'   in the release (`model_info()[[variant]]$sim_batchsize`). The stochastic
#'   draws depend on it, so keep it fixed when results are compared.
#' @param threads torch CPU threads; the default is the release's pinned value
#'   (1, the setting the models were fitted and validated with).
#' @param ... Further arguments for [FINN::simulateForest()].
#' @export
simulate_model = function(env, init_cohort = NULL, variant = NULL, batchsize = NULL, threads = NULL, ...) {
  info = model_info()
  if (is.null(variant)) variant = names(info)[1]
  v = info[[variant]]
  if (is.null(batchsize)) batchsize = v$sim_batchsize
  if (is.null(threads)) threads = v$threads
  if (!is.null(threads)) torch::torch_set_num_threads(as.integer(threads))
  args = utils::modifyList(list(patches = v$patches, patch_size = v$patch_size, batchsize = batchsize), list(...))
  do.call(FINN::simulateForest, c(list(load_model(variant), env = env, init_cohort = init_cohort), args))
}
