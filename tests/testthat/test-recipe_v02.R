## Round trips of build_model(): the zoo v0.2 recipe structure (R/recipe_v02.R)
## and the released v0.1 structure (R/definition.R).
##
## 1. Synthetic v0.2 round trip, always run: a recipe model with random weights
##    is dumped the way fit_stageB.R does and rebuilt through build_model();
##    every tensor must come back identical and a seeded simulateForest() of the
##    two models must agree exactly.
## 2. Released v0.2 round trip, run when FINN_ZOO_V02_FIT_DIR names a folder with
##    the recipe's parameters.pt and model.pt (a stage-B fit or smoke fit):
##    build_model(parameters.pt) against model.pt, tensor for tensor, plus one
##    seeded simulation.
## 3. Released v0.1 round trip, run when FINN_ZOO_V01_PARAMETERS names a v0.1
##    parameter file (and FINN_ZOO_V01_MODEL the matching model.pt, optional):
##    build_model(version = "v0.1") still loads every tensor and simulates.
## This file is the same in every FINN-zoo model package.

torch::torch_set_num_threads(1L)
quiet = function(expr) invisible(utils::capture.output(suppressMessages(expr)))

tensor_diff = function(a, b) {
  stopifnot(setequal(names(a), names(b)))
  vapply(names(a), function(nm) {
    x = a[[nm]]$detach()$cpu(); y = b[[nm]]$detach()$cpu()
    if (!identical(as.integer(x$shape), as.integer(y$shape))) return(Inf)
    if (x$numel() == 0) return(0)
    as.numeric((x$to(dtype = torch::torch_float64()) - y$to(dtype = torch::torch_float64()))$abs()$max())
  }, numeric(1))
}

## a small env table across the climate space: 2 sites x `years` rows, the climate
## at centre -/+ half a scale, the management column 0 / 1
test_env = function(env_scaling, env_vars, linear_var = NULL, years = 4L) {
  env_scaling = as.data.frame(env_scaling)
  env = expand.grid(year = seq_len(years), siteID = 1:2)[, c("siteID", "year")]
  for (v in env_vars) {
    c0 = env_scaling$center[match(v, env_scaling$variable)]; s0 = env_scaling$scale[match(v, env_scaling$variable)]
    env[[v]] = c0 + ifelse(env$siteID == 1L, -0.5, 0.5) * s0
  }
  if (!is.null(linear_var)) env[[linear_var]] = as.numeric(env$siteID == 2L)
  env
}

seeded_sim = function(model, env, patch_size, seed = 11L) {
  FINN::FINN.seed(seed)
  out = NULL
  quiet(out <- FINN::simulateForest(model, env = env, init_cohort = NULL, patches = 1L, patch_size = patch_size, device = "cpu"))
  as.data.frame(out$wide$site)
}

variants = if (is.null(model_info()$process)) names(model_info()) else "process"

for (variant in variants) {
  cfg = recipe_v02_config(variant)
  env_formula = stats::as.formula(paste("~", paste(cfg$env_vars, collapse = " + ")))

  test_that(paste0("v0.2 recipe: a synthetic model round-trips through build_model() [", variant, "]"), {
    n_species = 3L
    species = data.frame(species = seq_len(n_species), species_name = paste0("sp", seq_len(n_species)))
    env_scaling = data.frame(variable = cfg$env_vars, center = 0, scale = 1, stringsAsFactors = FALSE)
    env_scaling$center[match(cfg$env_vars, env_scaling$variable)] = seq_along(cfg$env_vars) * 10
    if (!is.null(cfg$linear_var)) { i = match(cfg$linear_var, env_scaling$variable); env_scaling$center[i] = 0; env_scaling$scale[i] = 1 }
    settings = list(prep_version = "synthetic", patch_size_ha = 0.05, n_patches = 1L, recruits_dbh = 10, reg_floor = 1e-3,
                    recruit_obs_weight = 0.05, env_formula = deparse(env_formula), linear_var = cfg$linear_var,
                    sen_max = cfg$sen_max, growth_scale = cfg$growth_scale, regeneration_sampler = cfg$regeneration_sampler,
                    theta = c(0.2, 5, 100)[seq_len(n_species)], recipe = "zoo v0.2 (synthetic test)")
    FINN::FINN.seed(3L)
    original = NULL
    quiet(original <- recipe_v02_structure(settings, n_species, env_scaling, config = cfg))
    ## the state the fit leaves behind: non-trivial values in every tensor the fit moves
    torch::with_no_grad({
      for (nm in c("reg_logK", "mort_sen_raw")) if (!is.null(original[[nm]])) original[[nm]]$copy_(torch::torch_randn(n_species) * 0.5)
      if (!is.null(original$nn_growth$beta)) original$nn_growth$beta$copy_(torch::torch_randn(n_species) * 0.1)
    })
    parameters = recipe_v02_parameters(original, species, settings)
    ## through torch_save / torch_load, as a release file travels
    f = tempfile(fileext = ".pt"); on.exit(unlink(f), add = TRUE)
    torch::torch_save(parameters, f)
    parameters = torch::torch_load(f)
    expect_true(is_recipe_v02(parameters))
    rebuilt = NULL
    quiet(rebuilt <- build_model(variant = variant, parameters = parameters, version = "v0.2-recipe"))
    quiet(auto <- build_model(variant = variant, parameters = parameters))
    d = tensor_diff(original$state_dict(), rebuilt$state_dict())
    expect_true(all(d == 0), info = paste("tensors differ:", paste(names(d)[d != 0], collapse = ", ")))
    expect_true(all(tensor_diff(original$state_dict(), auto$state_dict()) == 0))
    expect_identical(class(rebuilt$nn_growth)[1], if (is.null(cfg$linear_var)) "nn_sequential" else "env_dnn_linmgmt")
    expect_identical(formals(rebuilt$mortality_func)$growth_scale, cfg$growth_scale)
    expect_identical(rebuilt$light_mask, cfg$light_mask); expect_identical(rebuilt$n_quantiles, cfg$n_quantiles)
    expect_false(rebuilt$training)
    env = test_env(env_scaling, cfg$env_vars, cfg$linear_var)
    A = seeded_sim(original, env, settings$patch_size_ha); B = seeded_sim(rebuilt, env, settings$patch_size_ha)
    expect_equal(A, B, tolerance = 0)
    ## a second seeded run of the rebuilt model reproduces itself (the draws are seeded)
    expect_equal(B, seeded_sim(rebuilt, env, settings$patch_size_ha), tolerance = 0)
  })

  test_that(paste0("v0.2 recipe: build_model() refuses a parameter file that contradicts the yml [", variant, "]"), {
    n_species = 2L
    settings = list(patch_size_ha = 0.05, recruits_dbh = 10, reg_floor = 1e-3, env_formula = deparse(env_formula),
                    linear_var = cfg$linear_var, sen_max = cfg$sen_max, growth_scale = cfg$growth_scale, theta = c(1, 1),
                    recipe = "zoo v0.2 (synthetic test)")
    env_scaling = data.frame(variable = cfg$env_vars, center = 0, scale = 1, stringsAsFactors = FALSE)
    m = NULL
    quiet(m <- recipe_v02_structure(settings, n_species, env_scaling, config = cfg))
    p = recipe_v02_parameters(m, data.frame(species = 1:2), settings)
    p$state_dict[["par_growth_unconstrained"]] = NULL
    expect_error(quiet(build_model(variant = variant, parameters = p, version = "v0.2-recipe")), "missing in file: par_growth_unconstrained")
    p = recipe_v02_parameters(m, data.frame(species = 1:2), settings)
    p$settings$env_formula = "~ not_a_column"
    expect_error(quiet(build_model(variant = variant, parameters = p, version = "v0.2-recipe")), "env formula")
    p = recipe_v02_parameters(m, data.frame(species = 1:2), settings)
    p$settings$recipe = NULL
    expect_false(is_recipe_v02(p))
    expect_error(quiet(build_model(variant = variant, parameters = p, version = "v0.2-recipe")), "zoo v0.2")
  })
}

test_that("v0.2 recipe: build_model(parameters.pt) reproduces the recipe's model.pt (FINN_ZOO_V02_FIT_DIR)", {
  fit_dir = Sys.getenv("FINN_ZOO_V02_FIT_DIR", unset = "")
  skip_if(!nzchar(fit_dir), "FINN_ZOO_V02_FIT_DIR not set (a folder with the recipe's parameters.pt and model.pt)")
  variant = Sys.getenv("FINN_ZOO_V02_VARIANT", unset = names(model_info())[1])
  parameters = torch::torch_load(file.path(fit_dir, "parameters.pt"))
  expect_true(is_recipe_v02(parameters))
  saved = torch::torch_load(file.path(fit_dir, "model.pt"))
  rebuilt = NULL
  quiet(rebuilt <- build_model(variant = variant, parameters = parameters))
  d = tensor_diff(saved$state_dict(), rebuilt$state_dict())
  expect_true(all(d == 0), info = paste("tensors differ:", paste(names(d)[d != 0], collapse = ", ")))
  expect_identical(class(rebuilt$nn_growth)[1], class(saved$nn_growth)[1])
  expect_identical(rebuilt$light_mask, saved$light_mask); expect_identical(rebuilt$n_quantiles, saved$n_quantiles)
  expect_identical(rebuilt$regeneration_sampler, saved$regeneration_sampler); expect_identical(rebuilt$mortality_draw, saved$mortality_draw)
  expect_identical(formals(rebuilt$mortality_func)$growth_scale, formals(saved$mortality_func)$growth_scale)
  expect_equal(as.data.frame(rebuilt$env_scaling), as.data.frame(saved$env_scaling))
  cfg = recipe_v02_config(variant)
  env = test_env(parameters$env_scaling, cfg$env_vars, cfg$linear_var, years = 6L)
  ps = parameters$settings$patch_size_ha
  expect_equal(seeded_sim(saved, env, ps), seeded_sim(rebuilt, env, ps), tolerance = 0)
})

test_that("v0.1 structure: build_model(version = \"v0.1\") still loads a released parameter file (FINN_ZOO_V01_PARAMETERS)", {
  path = Sys.getenv("FINN_ZOO_V01_PARAMETERS", unset = "")
  skip_if(!nzchar(path), "FINN_ZOO_V01_PARAMETERS not set (a v0.1 parameters.pt)")
  variant = Sys.getenv("FINN_ZOO_V01_VARIANT", unset = names(model_info())[1])
  parameters = torch::torch_load(path)
  expect_false(is_recipe_v02(parameters))
  m = NULL
  quiet(m <- build_model(variant = variant, parameters = parameters, version = "v0.1"))
  quiet(m_auto <- build_model(variant = variant, parameters = parameters))
  expect_true(all(tensor_diff(m$state_dict(), m_auto$state_dict()) == 0))
  ## every tensor of the file that names a tensor of the model carries the model's value
  flat = if (!is.null(parameters$state_dict)) parameters$state_dict else parameters
  flat = flat[vapply(flat, inherits, logical(1), "torch_tensor")]
  names(flat) = sub("^buffer\\.", "", names(flat))
  sd = m$state_dict()
  common = intersect(names(flat), names(sd))
  expect_gte(length(common), 5L)
  d = tensor_diff(flat[common], sd[common])
  expect_true(all(d == 0), info = paste("tensors differ:", paste(names(d)[d != 0], collapse = ", ")))
  env_scaling = as.data.frame(m$env_scaling)
  env = test_env(env_scaling, env_scaling$variable, NULL, years = 4L)
  ps = model_info()[[variant]]$patch_size
  if (is.null(ps)) ps = 0.05
  ## a v0.1 definition may leave the env layers to FINN's lazy creation on the first
  ## simulation, which draws their initial weights from the seeded stream: warm up first
  quiet(FINN::simulateForest(m, env = env, init_cohort = NULL, patches = 1L, patch_size = ps, device = "cpu"))
  A = seeded_sim(m, env, ps)
  expect_equal(A, seeded_sim(m, env, ps), tolerance = 0)
  model_path = Sys.getenv("FINN_ZOO_V01_MODEL", unset = "")
  if (nzchar(model_path)) {
    saved = torch::torch_load(model_path)
    common = intersect(names(saved$state_dict()), names(sd))
    d = tensor_diff(saved$state_dict()[common], sd[common])
    expect_true(all(d == 0), info = paste("tensors differ from model.pt:", paste(names(d)[d != 0], collapse = ", ")))
    expect_equal(seeded_sim(saved, env, ps), A, tolerance = 0)
  }
})
