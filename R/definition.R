#' Build the model structure and set its fitted parameters
#'
#' The model definition of this package: the same `createProcess()` and
#' `finn()` calls as `calibration/fit.R`, followed by copying the fitted
#' parameters into the new object.
#'
#' @param variant Name of the model variant (see [model_info()]).
#' @param parameters Named list of parameter tensors, as saved by the
#'   calibration script.
#' @return A fitted `FINN` model.
#' @keywords internal
build_model = function(variant, parameters) {
  stop("build_model() is not written yet for this package.", call. = FALSE)
}
