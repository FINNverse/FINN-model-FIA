#' Definitions of the model variants
#'
#' One entry per variant: the species count, the climate covariates, the
#' process functions and the settings `build_model()` needs to recreate the
#' structure that `calibration/<variant>/02_fit.R` fitted. Kept in one place so
#' the quickstart, [simulate_model()] and the calibration scripts read the same
#' numbers.
#'
#' @return A named list of variant definitions.
#' @keywords internal
variant_definitions = function() {
  list(
    ## Oregon, 10-year time step, 27 species + other, linear climate response
    oregon = list(
      n_species    = 28L,
      covariates   = c("temp", "prec"),
      formula      = ~ temp + prec,
      hidden       = NULL,
      dropout      = 0,
      regeneration = "saturation",
      reg_K_init   = 500,
      adult_x_init = NULL,
      adult_A_init = NULL,
      reg_floor    = 0.01,
      recruits_dbh = 12.9,
      patches      = 4L,
      patch_size   = 0.01681,
      step_years   = 10L
    ),
    ## v0.2 Oregon (data oregon-10yr-v2): the species set S of 18 species under the
    ## site rule (species outside S < 2 % of basal area at the interval-start
    ## censuses), no lumped class; the recruitment cap K is fixed per species (not
    ## trained) and theta is frozen, both carried in the parameter list
    oregon_v02 = list(
      n_species    = 18L,
      covariates   = c("temp", "prec"),
      formula      = ~ temp + prec,
      hidden       = NULL,
      dropout      = 0,
      regeneration = "saturation",
      reg_K_init   = 500,
      adult_x_init = NULL,
      adult_A_init = NULL,
      reg_floor    = 0.01,
      recruits_dbh = 12.9,
      patches      = 4L,
      patch_size   = 0.01681,
      step_years   = 10L
    ),
    ## v0.2 conterminous US (data conus-annual-v2): the species set S of 36 species
    ## under the site rule (species outside S < 2 % of basal area at the interval-
    ## start censuses); the v0.1 structure otherwise; K fixed and theta frozen in
    ## the fit, both carried in the parameter list
    conus_v02 = list(
      n_species    = 36L,
      covariates   = c("tavg_c", "prec_mm", "cwd_z"),
      formula      = ~ tavg_c + prec_mm + cwd_z,
      hidden       = c(20L, 20L),
      dropout      = 0.1,
      regeneration = "adult",
      reg_K_init   = 50,
      adult_x_init = 20,
      adult_A_init = 15,
      reg_floor    = 1e-3,
      recruits_dbh = 12.9,
      patches      = 4L,
      patch_size   = 0.01681,
      step_years   = 1L
    ),
    ## conterminous US, annual time step, 20 species, climate networks
    conus = list(
      n_species    = 20L,
      covariates   = c("tavg_c", "prec_mm", "cwd_z"),
      formula      = ~ tavg_c + prec_mm + cwd_z,
      hidden       = c(20L, 20L),
      dropout      = 0.1,
      regeneration = "adult",
      reg_K_init   = 50,
      adult_x_init = 20,
      adult_A_init = 15,
      reg_floor    = 1e-3,
      recruits_dbh = 12.9,
      patches      = 4L,
      patch_size   = 0.01681,
      step_years   = 1L
    )
  )
}

#' Build the model structure and set its fitted parameters
#'
#' The model definition of this package: the same `createProcess()` and
#' `finn()` calls as `calibration/<variant>/02_fit.R`, followed by copying the
#' fitted parameters into the new object. FINN creates the climate networks
#' lazily on first use, so a one-site, one-step simulation on a dummy
#' environment is run before the parameters are copied.
#'
#' @param variant Name of the model variant (see [model_info()]).
#' @param parameters Named list of parameter tensors, as saved by the
#'   calibration script (`parameters.pt`).
#' @param env_scaling Data frame with `variable`, `center`, `scale`: the
#'   standardisation FINN learned from the training environment
#'   (`env_autoscale = TRUE`). `NULL` reads it from the package's
#'   `inst/model/covariates.csv`.
#' @return A fitted `FINN` model in evaluation mode.
#' @keywords internal
build_model = function(variant, parameters, env_scaling = NULL) {
  defs = variant_definitions()
  if (!variant %in% names(defs))
    stop("Unknown variant '", variant, "'; available: ", paste(names(defs), collapse = ", "), call. = FALSE)
  d = defs[[variant]]
  if (is.null(env_scaling)) {
    path = system.file("model", "covariates.csv", package = utils::packageName())
    if (!nzchar(path)) path = file.path("inst", "model", "covariates.csv")
    cov = utils::read.csv(path)
    cov = cov[cov$variant == variant, ]
    env_scaling = data.frame(variable = cov$name, center = cov$center, scale = cov$scale)
  }
  env_scaling = as.data.frame(env_scaling)[, c("variable", "center", "scale")]
  if (!setequal(env_scaling$variable, d$covariates))
    stop("env_scaling must cover exactly the covariates ", paste(d$covariates, collapse = ", "), call. = FALSE)

  ## createProcess() reads its formula from the call; through do.call() the
  ## formula object itself stands in the call, so the same code serves every variant
  process = function(func, custom_parameters = NULL)
    do.call(FINN::createProcess, c(list(formula = d$formula, func = func, hidden = d$hidden, dropout = d$dropout,
                                        optimizeSpecies = TRUE, optimizeEnv = TRUE),
                                   if (!is.null(custom_parameters)) list(custom_parameters = custom_parameters)))
  regeneration_process = switch(d$regeneration,
    saturation = process(FINN::regeneration_saturation,
                         list(reg_logK = rep(log(d$reg_K_init), d$n_species))),
    adult      = process(FINN::regeneration_adult,
                         list(reg_logK       = rep(log(d$reg_K_init),   d$n_species),
                              reg_adult_logx = rep(log(d$adult_x_init), d$n_species),
                              reg_adult_logA = rep(log(d$adult_A_init), d$n_species))),
    stop("unknown regeneration process ", d$regeneration))
  model = FINN::finn(
    N_species    = d$n_species,
    recruits_dbh = d$recruits_dbh,
    reg_floor    = d$reg_floor,
    competition_process  = FINN::createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
    growth_process       = process(FINN::growth),
    mortality_process    = process(FINN::mortality),
    regeneration_process = regeneration_process
  )
  model$env_scaling = env_scaling

  ## instantiate the lazily created climate networks
  dummy = data.frame(siteID = 1L, year = 1L)
  for (v in d$covariates) dummy[[v]] = env_scaling$center[match(v, env_scaling$variable)]
  invisible(FINN::simulateForest(model, env = dummy, init_cohort = NULL, patches = 1L,
                                 patch_size = d$patch_size, device = "cpu"))

  have = model$parameters
  missing = setdiff(names(have), names(parameters)); extra = setdiff(names(parameters), names(have))
  ## a value fixed in the fit (e.g. the v0.2 recruitment cap `reg_logK`, declared with
  ## trainable = FALSE) is a buffer of the model, not a parameter: it is carried in the
  ## list under its field name and copied into that buffer
  buffers = extra[vapply(extra, function(nm) inherits(model[[nm]], "torch_tensor"), logical(1))]
  extra = setdiff(extra, buffers)
  if (length(missing) || length(extra))
    stop("parameter list does not match the model: missing ", paste(missing, collapse = ", "),
         "; unexpected ", paste(extra, collapse = ", "), call. = FALSE)
  torch::with_no_grad({
    for (nm in names(have)) {
      if (!identical(as.integer(have[[nm]]$shape), as.integer(parameters[[nm]]$shape)))
        stop("parameter '", nm, "' has shape ", paste(parameters[[nm]]$shape, collapse = "x"),
             " but the model expects ", paste(have[[nm]]$shape, collapse = "x"), call. = FALSE)
      have[[nm]]$copy_(parameters[[nm]]$to(dtype = have[[nm]]$dtype))
    }
    for (nm in buffers) {
      if (!identical(as.integer(model[[nm]]$shape), as.integer(parameters[[nm]]$shape)))
        stop("buffer '", nm, "' has shape ", paste(parameters[[nm]]$shape, collapse = "x"),
             " but the model expects ", paste(model[[nm]]$shape, collapse = "x"), call. = FALSE)
      model[[nm]]$copy_(parameters[[nm]]$to(dtype = model[[nm]]$dtype))
    }
  })
  model$eval()
  model
}
