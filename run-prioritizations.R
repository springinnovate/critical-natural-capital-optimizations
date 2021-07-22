library(sf)
library(raster)
library(Matrix)
library(tidyverse)
library(prioritizr)
library(gurobi)
library(doParallel)
library(foreach)
select <- dplyr::select
message_time <- function(...) {
  submessage <- paste(...)
  message(glue::glue("({Sys.time()}) {submessage}"))
}

n_cores <- 10
parallel_threshold_km2 <- 3e6

data_dir <- "data/"
rij_dir <- file.path(data_dir, "prioritizr", "rij_2km")
out_dir <- "solutions/"

# features
features <- file.path(data_dir, "prioritizr", "features_2km.csv") %>% 
  read_csv()

# run matrix
runs <- file.path(data_dir, "optimization-runs.csv") %>% 
  read_csv() %>% 
  nest(layers = layer)

# raster template
r_pu <- file.path(data_dir, "pu", "pu_eck4_2km.tif") %>% 
  raster() %>% 
  raster()

# maximum planning unit number
pu_lookup <- str_glue("pu_2km.rds") %>% 
  file.path("data", "prioritizr", .) %>% 
  readRDS()
n_pu <- nrow(pu_lookup)

# run each prioritization scenario in sequence
solution_df <- NULL
for (scn in seq_len(nrow(runs))) {
  run_scn <- runs[scn, ]
  scn_dir <- file.path(out_dir, run_scn$scenario)
  dir.create(scn_dir, recursive = TRUE, showWarnings = FALSE)
  
  message_time("Scenario:", run_scn$scenario)
  
  # targets to use for this scenario
  if (isTRUE(run_scn$variable_target)) {
    targets <- seq(5, 100, by = 5)
  } else {
    targets <- 90
  }
  
  # expand to all features, global or country-level
  feature_scn <- unnest(run_scn, cols = layers)
  if (isTRUE(run_scn$country_prioritization)) {
    feature_scn <- features %>% 
      filter(region != "global") %>% 
      select(layer, region, name) %>% 
      inner_join(feature_scn, ., by = "layer")
  } else {
    feature_scn <- features %>% 
      filter(region == "global") %>% 
      select(layer, region, name) %>% 
      inner_join(feature_scn, ., by = "layer") %>% 
      # fix name for global rij matrices
      mutate(name = str_replace(name, "^global_", "000_"))
  }
  regions <- unique(feature_scn$region) %>% 
    sort()
  
  # prioritize regions independently
  for (rgn in regions) {
    message_time("Region:", rgn)
    rgn_code <- if_else(rgn == "global", "000", rgn)
    
    # load features
    feature_rgn <- feature_scn %>% 
      filter(region == rgn) %>% 
      mutate(id = row_number()) %>% 
      select(id, name)
    rij <- feature_rgn %>% 
      mutate(rij = str_glue("rij_{name}.rds") %>% 
               file.path(rij_dir, .) %>% 
               map(read_rds)) %>% 
      select(id, rij) %>% 
      unnest(rij)
    
    # create sparse matrix
    rij_mat <- sparseMatrix(i = rij$id, j = rij$pu, x = rij$amount,
                            index1 = TRUE, use.last.ij = FALSE,
                            dims = c(max(rij$id), n_pu))
    # remove planning units not within region
    pu_included <- sort(unique(rij$pu))
    # connects planning units for this region to cell ids in the global raster
    pu_cell_id_rgn <- pu_lookup$cell_id[pu_included]
    rij_mat <- rij_mat[, pu_included, drop = FALSE]
    n_pu_rgn <- ncol(rij_mat)
    message_time("Priotizing", nrow(rij_mat), "features in", n_pu_rgn, 
                 "planning units")
    
    # construct problem, use cost = 1
    # leave target missing for now
    p <- problem(rep(1, n_pu_rgn), 
                 features = feature_rgn, 
                 rij_matrix = rij_mat,
                 run_checks = FALSE) %>%
      add_min_set_objective() %>%
      add_binary_decisions()
    rm(rij_mat, rij)
    
    # paralellize over targets for small countries
    rgn_area_km2 <- prod(res(r_pu) / 1000) * n_pu_rgn
    if (rgn_area_km2 > parallel_threshold_km2) {
      registerDoParallel(cores = 1)
      rgn_cores <- n_cores
    } else {
      registerDoParallel(cores = n_cores)
      rgn_cores <- 1
    }
    # apply different targets
    sol_targets <- foreach (t = targets, .combine = bind_rows) %dopar% {
      message_time("Prioritizing: Scenario", run_scn$scenario, "- Region", rgn, 
                   "- Target", t)
      # set target
      p <- add_relative_targets(p, t / 100) %>% 
        add_gurobi_solver(gap = 0.01, threads = rgn_cores)
      
      # solve
      # need to turn off checks due to large range of features
      s <- solve(p, run_checks = FALSE)
      
      # convert to raster
      r_sol <- r_pu
      r_sol[pu_cell_id_rgn] <- as.integer(s)
      f_out <- str_glue("solution_scenario-{run_scn$scenario}_",
                        "{rgn_code}_target-{t}.tif") %>% 
        file.path(scn_dir, .)
      r_sol <- writeRaster(r_sol, filename = f_out, overwrite = TRUE,
                           options = c("COMPRESS=DEFLATE"))
      
      # compile solution data
      tibble(
        scenario = run_scn$scenario,
        region = rgn,
        target = t,
        objective = attr(s, "objective"),
        prop_selected_available = attr(s, "objective") / n_pu_rgn,
        prop_selected_total = attr(s, "objective") / n_pu,
        runtime = attr(s, "runtime"),
        solution_raster = basename(f_out)
      )
    }
    solution_df <- bind_rows(solution_df, sol_targets)
  }
}
write_csv(solution_df, file.path(out_dir, "solutions.csv"))