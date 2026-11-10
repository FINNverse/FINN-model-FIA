## The zoo v0.2 recipe: the released model structure, rebuilt from parameters.pt
##
## This file is the same in every FINN-zoo model package (the recipe is
## FINN-zoo-data/RECIPE_v02.md). It holds the recipe's custom mortality kernel
## and environment module (sections 3.1 and 3.2 of the recipe, verbatim) and
## build_model_v02(), which recreates the stage-B simulation model from the
## released `parameters.pt` (state_dict + env_scaling + species + settings) and
## the package's `inst/model/model_v02.yml` (the recipe constants the parameter
## file does not carry). Change it in one package, then copy it to the others;
## what differs between the packages is the yml, not this file.

utils::globalVariables("self")

#' U-shaped mortality of the zoo v0.2 recipe
#'
#' `FINN::mortality`'s linear predictor (light, size, growth in percent) plus a
#' size term `sen * (dbh / 100)^2` with `sen = sen_max * sigmoid(mort_sen_raw)`
#' per species in `[0, sen_max]`. The two custom parameters are registered on
#' the model by `createProcess(custom_parameters = )`: `mort_sen_raw` (fitted
#' in stage A, frozen afterwards) and `mort_sen_max` (fixed). The `growth_scale`
#' formal is the one `FINN::mortality` has, so FINN initialises the growth
#' coefficient on the same scale; the released value is set through
#' `createProcess(growth_scale = )`.
#'
#' @inheritParams FINN::mortality
#' @keywords internal
mortality_ushape = function(dbh, species, trees, parMort, pred, light, base_steepness = 5, debug = FALSE, growth = NULL, growth_scale = 0.01) {
  if (is.null(growth)) growth <- self$g
  d   <- dbh / 100
  sen <- (self$mort_sen_max * torch::torch_sigmoid(self$mort_sen_raw))[species]
  torch::torch_sigmoid(pred + growth / growth_scale * parMort[, 3][species] + light * parMort[, 1][species] +
                       parMort[, 2][species] * d + sen * d * d)
}

#' Environment module of the zoo v0.2 recipe: climate network + linear management
#'
#' `env(x) = MLP(climate) + beta_s * mgmt`. FINN's env design matrix is
#' `stats::model.matrix(formula, env)`: the intercept column first, then the
#' formula's terms in order; `mgmt_col` is the management column's position in
#' it. `mlp` is FINN's own network layout (`build_NN`: linear, SELU, dropout
#' after every hidden activation, linear output with bias) on every other column
#' (the constant intercept column stays in, as in FINN's networks); `beta` is
#' one linear weight per species on the management column, added outside the
#' network, so the flag keeps its own scale and `beta_s` is the shift of a
#' managed plot on the log-growth / mortality logit. The column indices are
#' plain R fields (kept by `torch_save()` with the module). Passed to
#' `createProcess(NN = )`; the state_dict names are `nn_<process>.mlp.<i>.*`
#' and `nn_<process>.beta`.
#'
#' @keywords internal
env_dnn_linmgmt = torch::nn_module("env_dnn_linmgmt",
  initialize = function(n_inputs, n_species, mgmt_col, hidden, dropout, beta_init) {
    self$mgmt_col  <- as.integer(mgmt_col)
    self$clim_cols <- setdiff(seq_len(n_inputs), self$mgmt_col)
    layers <- list(); n_in <- length(self$clim_cols)
    for (h in hidden) {
      layers <- c(layers, torch::nn_linear(n_in, h), torch::nn_selu())
      if (dropout > 0) layers <- c(layers, torch::nn_dropout(p = dropout))
      n_in <- h
    }
    self$mlp  <- do.call(torch::nn_sequential, c(layers, torch::nn_linear(n_in, n_species)))
    self$beta <- torch::nn_parameter(torch::torch_tensor(beta_init, dtype = torch::torch_float32()))
  },
  forward = function(x) {
    self$mlp(x[.., self$clim_cols]) + x[.., self$mgmt_col, drop = FALSE] * self$beta
  })

#' Is this parameter list a zoo v0.2 recipe model?
#'
#' The recipe's `fit_stageB.R` tags `parameters$settings$recipe` with a string
#' that starts with `"zoo v0.2"`; a v0.1 parameter file (a flat list of tensors,
#' or a list whose settings carry no such tag) is not.
#'
#' @param parameters The object `torch::torch_load()` returns for a parameter file.
#' @return `TRUE` or `FALSE`.
#' @keywords internal
is_recipe_v02 = function(parameters) {
  is.list(parameters) && !is.null(parameters$settings) && is.list(parameters$settings) &&
    is.character(parameters$settings$recipe) && grepl("^zoo v0\\.2", parameters$settings$recipe[1])
}

#' The recipe constants of this package: `inst/model/model_v02.yml`
#'
#' What `parameters.pt` does not carry: the env columns and the management
#' column (checked against the saved formula), the network layout (checked
#' against the state_dict shapes), the competition settings, the mortality draw,
#' the recruitment sampler and the regeneration process. Entries of the yml's
#' `variants:` block override the top level for that variant (FIA: oregon and
#' conus differ in their env columns).
#'
#' @param variant Model variant; `NULL` takes the top-level entries only.
#' @return A named list.
#' @keywords internal
recipe_v02_config = function(variant = NULL) {
  path = system.file("model", "model_v02.yml", package = utils::packageName())
  if (!nzchar(path)) stop("inst/model/model_v02.yml is missing from this package.", call. = FALSE)
  cfg = yaml::read_yaml(path)
  if (!is.null(variant) && !is.null(cfg$variants)) {
    if (!is.null(cfg$variants[[variant]])) cfg = utils::modifyList(cfg, cfg$variants[[variant]])
  }
  cfg$variants = NULL
  if (!identical(cfg$structure, "v0.2-recipe"))
    stop("model_v02.yml: `structure` must be \"v0.2-recipe\" (found ", deparse(cfg$structure), ").", call. = FALSE)
  for (nm in c("env_vars", "hidden", "dropout", "light_mask", "n_quantiles", "mortality_draw", "regeneration_sampler",
               "regeneration", "sen_max", "growth_scale"))
    if (is.null(cfg[[nm]]) && nm != "hidden") stop("model_v02.yml lacks `", nm, "`.", call. = FALSE)
  cfg$env_vars    = as.character(unlist(cfg$env_vars))
  cfg$hidden      = as.integer(unlist(cfg$hidden))          # NULL or [] = a linear env layer
  cfg$n_quantiles = as.integer(cfg$n_quantiles)
  cfg$light_mask  = isTRUE(cfg$light_mask)
  if (!is.null(cfg$linear_var) && !nzchar(cfg$linear_var)) cfg$linear_var = NULL
  if (!cfg$regeneration %in% c("saturation", "frozen_zero"))
    stop("model_v02.yml: `regeneration` must be \"saturation\" or \"frozen_zero\".", call. = FALSE)
  cfg
}

#' The network layout a recipe state_dict describes
#'
#' Reads, for one process, whether the env response is the recipe's management
#' module (`nn_<type>.mlp.<i>.*` + `nn_<type>.beta`) or FINN's own network
#' (`nn_<type>.<i>.*`), the hidden layer widths (the output size of every
#' linear layer but the last) and whether dropout layers sit between them
#' (FINN's `build_NN` and the module place linear, SELU, dropout, so the linear
#' layers are 3 apart with dropout and 2 apart without).
#'
#' @param state_dict Named list of tensors.
#' @param type `"growth"`, `"mortality"` or `"regeneration"`.
#' @return A list: `module`, `hidden`, `dropout_layers`, `n_inputs`.
#' @keywords internal
recipe_v02_layout = function(state_dict, type) {
  prefix = paste0("nn_", type, ".")
  nms = grep(paste0("^", prefix), names(state_dict), value = TRUE)
  if (!length(nms)) stop("The state_dict holds no tensors of nn_", type, ".", call. = FALSE)
  module = any(startsWith(nms, paste0(prefix, "mlp.")))
  core = sub(paste0("^", prefix, if (module) "mlp\\." else ""), "", nms)
  w = core[grepl("^[0-9]+\\.weight$", core)]
  idx = sort(as.integer(sub("\\.weight$", "", w)))
  if (!length(idx)) stop("The state_dict holds no linear layer of nn_", type, ".", call. = FALSE)
  shape_of = function(i) as.integer(state_dict[[paste0(prefix, if (module) "mlp." else "", i, ".weight")]]$shape)
  list(module = module,
       hidden = if (length(idx) > 1) vapply(idx[-length(idx)], function(i) shape_of(i)[1], integer(1)) else integer(0),
       dropout_layers = length(idx) > 1 && all(diff(idx) == 3L),
       n_inputs = shape_of(idx[1])[2])
}

#' The recipe model structure, before the parameters are set
#'
#' The `createProcess()` and `finn()` calls of the recipe's `fit_stageB.R`
#' (competition, growth and mortality frozen; `regeneration_saturation` with the
#' frozen dispersion, or the frozen zero recruitment of the prescribed variant),
#' followed by the dummy simulation that makes FINN create the env networks.
#' Also used by the tests to make a synthetic recipe model.
#'
#' @param settings `parameters$settings` of the recipe (at least `env_formula`,
#'   `recruits_dbh`, `reg_floor`, `patch_size_ha`; `linear_var`, `sen_max`,
#'   `growth_scale`, `regeneration_sampler`, `theta`, `recruit_obs_weight` with
#'   the yml / FINN defaults as fallback).
#' @param n_species Number of species.
#' @param env_scaling Data frame `variable`, `center`, `scale` covering the env
#'   columns (the management column with centre 0 and scale 1).
#' @param config [recipe_v02_config()].
#' @param layout Network layout: `hidden` (integer vector), `dropout_layers`
#'   (logical), `module` (logical: the management module); the default is the
#'   yml's.
#' @return A `FINN` model with randomly initialised parameters and its env
#'   networks created, in eval mode.
#' @keywords internal
recipe_v02_structure = function(settings, n_species, env_scaling, config = recipe_v02_config(),
                                layout = list(hidden = config$hidden, dropout_layers = config$dropout > 0,
                                              module = !is.null(config$linear_var))) {
  s = settings
  pick = function(nm, default) if (!is.null(s[[nm]])) s[[nm]] else default
  env_formula = stats::as.formula(paste(s$env_formula, collapse = " "))
  env_vars = all.vars(env_formula)
  if (!identical(env_vars, config$env_vars))
    stop("The saved parameters were fitted with env formula `", paste(deparse(env_formula), collapse = ""),
         "`; model_v02.yml lists the env columns ", paste(config$env_vars, collapse = ", "), ".", call. = FALSE)
  linear_var = pick("linear_var", config$linear_var)
  if (!is.null(linear_var) && !linear_var %in% env_vars)
    stop("The management column `", linear_var, "` is not in the env formula.", call. = FALSE)
  if (isTRUE(layout$module) != !is.null(linear_var))
    stop(if (isTRUE(layout$module)) "The state_dict holds the management module (nn_growth.mlp.*, nn_growth.beta) but no management column is declared (settings$linear_var, model_v02.yml linear_var)."
         else "A management column is declared but the state_dict holds FINN's plain networks, not the env_dnn_linmgmt module.", call. = FALSE)
  env_scaling = as.data.frame(env_scaling, stringsAsFactors = FALSE)
  if (!all(env_vars %in% env_scaling$variable))
    stop("env_scaling lacks ", paste(setdiff(env_vars, env_scaling$variable), collapse = ", "), ".", call. = FALSE)
  hidden = as.integer(layout$hidden)
  if (!identical(hidden, config$hidden))
    stop("The state_dict describes hidden layers c(", paste(hidden, collapse = ", "), "); model_v02.yml says c(",
         paste(config$hidden, collapse = ", "), ").", call. = FALSE)
  ## the dropout layers must exist for the state_dict indices to match (linear 0, 3, 6);
  ## their rate is set to 0 after loading, as in the released (frozen) model
  dropout_build = if (isTRUE(layout$dropout_layers)) max(config$dropout, 1e-6) else 0
  sen_max = pick("sen_max", config$sen_max)
  growth_scale = pick("growth_scale", config$growth_scale)
  theta = pick("theta", rep(1, n_species))
  if (length(theta) == 1L) theta = rep(theta, n_species)

  ## the dummy env row: every climate column at its centre, management 0; also the
  ## design-matrix column names the module reads its management column from
  dummy = data.frame(siteID = 1L, year = 1L)
  for (v in env_vars) dummy[[v]] = env_scaling$center[match(v, env_scaling$variable)]
  if (!is.null(linear_var)) dummy[[linear_var]] = 0
  X_names = colnames(stats::model.matrix(env_formula, dummy[, env_vars, drop = FALSE]))
  new_module = function() env_dnn_linmgmt(n_inputs = length(X_names), n_species = n_species, mgmt_col = match(linear_var, X_names),
                                          hidden = hidden, dropout = dropout_build, beta_init = rep(0, n_species))
  ## createProcess() reads its formula from the call; through do.call() the formula
  ## object itself stands in the call, so one code path serves every package
  process = function(func, ...) {
    args = c(list(formula = env_formula, func = func), list(...))
    if (!is.null(linear_var)) args$NN = new_module()
    else { args$hidden = if (length(hidden)) hidden else NULL; args$dropout = dropout_build }
    do.call(FINN::createProcess, args)
  }
  regeneration_process = switch(config$regeneration,
    saturation = do.call(FINN::createProcess, list(
      formula = env_formula, func = FINN::regeneration_saturation, dispersion_parameter = theta,
      custom_parameters = list(reg_logK = list(init = rep(0, n_species), trainable = TRUE)),
      optimizeSpecies = TRUE, optimizeEnv = TRUE, regeneration_sampler = pick("regeneration_sampler", config$regeneration_sampler))),
    frozen_zero = do.call(FINN::createProcess, list(
      formula = env_formula, func = FINN::regeneration, hidden = NULL,
      optimizeSpecies = FALSE, optimizeEnv = FALSE, sample_regeneration = FALSE)))
  model = FINN::finn(
    N_species    = n_species,
    recruits_dbh = s$recruits_dbh,
    reg_floor    = s$reg_floor,
    recruit_obs_weight = pick("recruit_obs_weight", 1),
    competition_process  = FINN::createProcess(~0, FINN::competition, optimizeSpecies = FALSE,
                                               n_quantiles = config$n_quantiles, light_mask = config$light_mask),
    growth_process       = process(FINN::growth, optimizeSpecies = FALSE, optimizeEnv = FALSE),
    mortality_process    = process(mortality_ushape, optimizeSpecies = FALSE, optimizeEnv = FALSE,
                                   growth_scale = growth_scale, mortality_draw = config$mortality_draw,
                                   custom_parameters = list(mort_sen_raw = list(init = rep(0, n_species), trainable = FALSE),
                                                            mort_sen_max = list(init = sen_max, trainable = FALSE))),
    regeneration_process = regeneration_process
  )
  model$period_mode = "per_site"
  model$env_scaling = env_scaling
  ## FINN creates the env networks on the first simulation; a one-site, one-year
  ## bare-ground run creates them (and makes the modules registered submodules)
  invisible(utils::capture.output(suppressMessages(
    FINN::simulateForest(model, env = dummy, patches = 1L, patch_size = s$patch_size_ha, init_cohort = NULL, device = "cpu"))))
  if (!identical(formals(model$mortality_func)$growth_scale, growth_scale))
    stop("growth_scale was not set on the mortality kernel.", call. = FALSE)
  recipe_v02_freeze(model)
  model
}

#' Put a recipe model into its released state
#'
#' The frozen networks run without dropout and outside the gradient, the
#' dispersion too; eval mode.
#' @param model A `FINN` model built by [recipe_v02_structure()].
#' @keywords internal
recipe_v02_freeze = function(model) {
  for (nn in list(model$nn_growth, model$nn_mortality)) {
    for (m in nn$modules) if (inherits(m, "nn_dropout")) m$p = 0
    for (p in nn$parameters) invisible(p$requires_grad_(FALSE))
  }
  if (!is.null(model$par_theta_recruits_raw)) invisible(model$par_theta_recruits_raw$requires_grad_(FALSE))
  model$eval()
  invisible(model)
}

#' Build a zoo v0.2 recipe model from its released parameter file
#'
#' Rebuilds the stage-B simulation model of the recipe (`calibration/v02/recipe/
#' fit_stageB.R`): competition, `FINN::growth` and [mortality_ushape()] with the
#' climate networks (or the [env_dnn_linmgmt] modules of a model with a
#' management column), `FINN::regeneration_saturation` with the frozen
#' dispersion and the fitted cap `reg_logK`; then loads every tensor of
#' `parameters$state_dict` by name. The structure is read from the parameter
#' file itself (the env formula, the management column, the network layout from
#' the tensor names and shapes, theta, `sen_max`, `growth_scale`, the sampler)
#' and the rest from `inst/model/model_v02.yml`; any disagreement between the
#' two, any tensor missing, unexpected or of another shape is an error, never a
#' silent partial load.
#'
#' @param parameters Named list as saved by the recipe's `fit_stageB.R`:
#'   `state_dict` (every parameter and buffer tensor of the model, by name),
#'   `env_scaling` (centre and scale of the env columns), `species`, `settings`.
#' @param variant Model variant, for the yml's `variants:` block; `NULL` for a
#'   package with one structure.
#' @return A fitted `FINN` model in eval mode, ready for [FINN::simulateForest()].
#' @keywords internal
build_model_v02 = function(parameters, variant = NULL) {
  for (nm in c("state_dict", "env_scaling", "species", "settings"))
    if (is.null(parameters[[nm]])) stop("parameters.pt lacks `", nm, "`; not a zoo v0.2 recipe parameter file.", call. = FALSE)
  if (!is_recipe_v02(parameters))
    stop("parameters$settings$recipe is not a zoo v0.2 tag: ", deparse(parameters$settings$recipe), call. = FALSE)
  config = recipe_v02_config(variant)
  species = as.data.frame(parameters$species, stringsAsFactors = FALSE)
  n_species = nrow(species)
  want = parameters$state_dict
  lay_g = recipe_v02_layout(want, "growth"); lay_m = recipe_v02_layout(want, "mortality")
  if (!identical(lay_g$hidden, lay_m$hidden) || !identical(lay_g$module, lay_m$module) || !identical(lay_g$dropout_layers, lay_m$dropout_layers))
    stop("The growth and mortality networks of the state_dict differ in layout; not a recipe model.", call. = FALSE)
  model = recipe_v02_structure(parameters$settings, n_species, as.data.frame(parameters$env_scaling, stringsAsFactors = FALSE),
                               config = config, layout = lay_g)
  ## match the saved tensors to the rebuilt model by name; any difference in the
  ## set of names or in a tensor's shape is an error
  have = model$state_dict()
  if (!setequal(names(have), names(want)))
    stop("Parameter names of the saved file and the rebuilt model differ.\n  missing in file: ",
         paste(setdiff(names(have), names(want)), collapse = ", "), "\n  unknown in file: ",
         paste(setdiff(names(want), names(have)), collapse = ", "), call. = FALSE)
  bad = names(want)[!vapply(names(want), function(n) identical(as.integer(want[[n]]$shape), as.integer(have[[n]]$shape)), logical(1))]
  if (length(bad)) stop("Tensor shapes differ for: ", paste(bad, collapse = ", "), call. = FALSE)
  model$load_state_dict(want)
  recipe_v02_freeze(model)
  model
}

#' Dump a recipe model the way `fit_stageB.R` does
#'
#' The parameter list of a release: `state_dict`, `env_scaling`, `species`,
#' `settings`. Used by the tests; the calibration scripts write the same list.
#' @param model A recipe model.
#' @param species Data frame with a `species` column (and the names).
#' @param settings The settings list (with the `recipe` tag).
#' @keywords internal
recipe_v02_parameters = function(model, species, settings) {
  list(state_dict = model$state_dict(), env_scaling = model$env_scaling, species = species, settings = settings)
}
