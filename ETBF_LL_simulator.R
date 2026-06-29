# =============================================================================
# ETBF longline data simulator - DRIVER script.
#
# Runs simdatfun() (see simdatfun.R) to generate synthetic longline catch-and-
# effort data sets with a KNOWN true biomass trajectory, for testing CPUE
# standardisation / spatial models. This script: sets parameters, defines the
# catchability matrix qq, runs several scenarios, writes the data to CSV, and
# draws a diagnostic comparing nominal CPUE against the true biomass.
#
# Original: Simon Hoyle, Jan 2016 (after H. Winker). Revised for portability
# (no Windows-only calls, no hardcoded absolute paths) and reproducibility.
#
# To run: edit base_dir and fun_path below, then source this file.
# =============================================================================

## --- Configuration ------------------------------------------------------
# All output is written under base_dir; fun_path points to the function file.
base_dir <- file.path(getwd(), "sim_output")    # output root (change as needed)
fun_path <- "simdatfun.R"                         # location of simdatfun()
set.seed(42)                                      # makes every run reproducible

source(fun_path)        # defines simdatfun(); also loads mgcv, VGAM, boot
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

## --- Run constants ------------------------------------------------------
# These are used by the calls below and by the helper functions (via scoping).
nclus    <- 4           # unused by simdatfun(), kept for call compatibility
nspec    <- 10          # number of species
nyr      <- 20          # number of years
ntactics <- 4           # number of fishing tactics (rows of qq)
p        <- 1.3         # Tweedie power parameter
disp     <- 10          # dispersion (Tweedie phi / NB size)
spnames  <- c("ALB","BET","YFT","SWO","MLS","DOL","SBT","BUM","BLM","SHK")

# qq[tactic, species] = catchability of each species under each tactic.
# Tactics: 1 deep longline (BET), 2 shallow longline (SWO), 3 YFT+ALB, 4 other.
# Must have ntactics rows and nspec columns.
qq <- rbind(
  c(0.50,1.00,0.80,0.05,0.01,0.01,0.01,0.01,0.01,0.30),
  c(0.20,0.10,0.50,1.00,0.50,0.01,0.01,0.60,0.60,0.60),
  c(1.00,0.40,1.00,0.10,1.00,0.05,0.20,0.20,0.30,0.30),
  c(0.70,0.50,0.01,0.10,0.10,0.50,1.00,0.20,0.20,0.60))
colnames(qq) <- spnames

## --- Example single runs ------------------------------------------------
# A few illustrative data sets that vary grid size, dispersion, spatial model
# (spdist) and observation distribution (rtype). Each returns list(sets, B).
sets1 <- simdatfun(nclus, nspec, spnames, nyr=20, p=1.3, disp=10, qq=qq,
                   nvess=16, nlat=10, nlon=10, ntactics=4, tripmode=10,
                   spdist="ew_rndom", rtype="Tweedie")
sets2 <- simdatfun(nclus, nspec, spnames, nyr=20, p=1.3, disp=5,  qq=qq,
                   nvess=4,  nlat=4,  nlon=4,  ntactics=4, tripmode=3,
                   spdist="ew_trend", rtype="zinb")
sets3 <- simdatfun(nclus, nspec, spnames, nyr=10, p=1.3, disp=5,  qq=qq,
                   nvess=3,  nlat=4,  nlon=4,  ntactics=4, tripmode=3,
                   spdist="ew_trend", rtype="Tweedie")
# $sets = one row per set (effort + catch); $B = the true biomass matrix.
write.csv(sets2$sets, file.path(base_dir, "smallsim.csv"))
write.csv(sets2$B,    file.path(base_dir, "smallsimB.csv"))

## --- Replicate datasets (two spatial scenarios) -------------------------
# Write nrep independent replicates for a given spatial model. Uses the run
# constants and qq defined above. simData<i>.csv = sets, simB<i>.csv = biomass.
nrep <- 100
run_replicates <- function(out_dir, spdist) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (i in 1:nrep) {
    s <- simdatfun(nclus, nspec, spnames, nyr=20, p=1.3, disp=10, qq=qq,
                   nvess=16, nlat=10, nlon=10, ntactics=4, tripmode=10,
                   spdist=spdist, rtype="Tweedie")
    write.csv(s$sets, file.path(out_dir, paste0("simData", i, ".csv")))
    write.csv(s$B,    file.path(out_dir, paste0("simB", i, ".csv")))
  }
}
run_replicates(file.path(base_dir, "SDHsim3"),  "ew_trend")
run_replicates(file.path(base_dir, "SDHsim4c"), "ew_rndom")

## --- Scaling / memory-test datasets -------------------------------------
# Increasingly large data sets (names approximate the row count) for checking
# performance and memory use. nyr/nvess grow; other settings are held fixed.
mem_dir <- file.path(base_dir, "memtests")
dir.create(mem_dir, recursive = TRUE, showWarnings = FALSE)
mem_runs <- list(sim2000  = list(nyr=18, nvess=4),
                 sim4000  = list(nyr=36, nvess=4),
                 sim8000  = list(nyr=36, nvess=8),
                 sim16000 = list(nyr=36, nvess=16))
for (nm in names(mem_runs)) {
  cfg <- mem_runs[[nm]]
  s <- simdatfun(nclus, nspec, spnames, nyr=cfg$nyr, p=1.3, disp=5, qq=qq,
                 nvess=cfg$nvess, nlat=6, nlon=6, ntactics=4, tripmode=3,
                 spdist="ew_trend", rtype="Tweedie")
  write.csv(s$sets, file.path(mem_dir, paste0(nm, ".csv")))
  write.csv(s$B,    file.path(mem_dir, paste0(nm, "B.csv")))
}

## --- Diagnostic: nominal CPUE vs true biomass ---------------------------
# For each species, plot the true biomass against nominal CPUE - the year-mean
# catch, either expected ("_M") or observed ("_R") - with both normalised to
# their own mean so the trends are comparable. If standardisation worked
# perfectly the dashed line would track the points. Writes a PNG (portable;
# replaces the original Windows-only windows()/savePlot()). Adapts to any nspec.
plot_cpue_vs_B <- function(res, type = "M",
                           file = file.path(base_dir, paste0("CPUE_vs_B_", type, ".png"))) {
  dat <- res$sets; B <- res$B
  nyr <- nrow(B); nsp <- ncol(B)
  cols <- grep(paste0("_", type, "$"), names(dat))          # the _M (or _R) columns
  spnm <- sub(paste0("_", type, "$"), "", names(dat)[cols]) # species names
  stdz <- function(j) { a <- tapply(dat[, j], dat$yr, mean); a / mean(a) }  # year-mean / mean
  nc <- ceiling(sqrt(nsp)); nr <- ceiling(nsp / nc)         # panel grid
  png(file, width = 1200, height = 900); on.exit(dev.off()) # close device on exit
  par(mfrow = c(nr, nc), mar = c(3, 3, 2, 1))
  for (i in seq_len(nsp)) {
    plot(1:nyr, B[, i] / mean(B[, i]), ylim = c(0, 2.5), pch = 1,
         main = spnm[i], xlab = "year", ylab = "rel. index")  # points = true biomass
    lines(1:nyr, stdz(cols[i]), lty = 2)                      # dashed = nominal CPUE
  }
}
plot_cpue_vs_B(sets1, "M")   # expected catch vs biomass
plot_cpue_vs_B(sets1, "R")   # observed catch vs biomass
