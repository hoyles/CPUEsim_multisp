# =============================================================================
# simdatfun() - simulate longline catch-and-effort data for the Eastern Tuna
# and Billfish Fishery (ETBF), with a KNOWN "true" biomass trajectory.
#
# Why: produce synthetic data whose underlying abundance is known exactly, so
# CPUE-standardisation and spatial models can be tested against the truth.
#
# Original: Simon Hoyle, Jan 2016 (after H. Winker).
# Revision: bug fixes (tactic join, nsets alignment, month guard, do_ewtrend
#   axis, list-vs-matrix robustness) + generalised to any nspec / ntactics +
#   vessel-q log-mean 0 + density normalised to mean 1.
#
# CATCH MODEL (per longline set):
#   expected catch(species) = effort x q x density x biomass
#     effort  = 1                                  (each set = one unit of effort)
#     q       = q_tactic x q_vessel x q_month       (catchability)
#     density = the species' relative density in the set's grid cell
#     biomass = the species' true biomass that year
#   Observed catch is then drawn from a Tweedie or zero-inflated NB distribution.
#
# SAMPLING HIERARCHY (each level nested in the one above):
#   vessel  ->  vessel-year  ->  trip  ->  set       (set = unit of observation)
#
# DEPENDENCIES: mgcv (rTweedie), VGAM (rzinegbin), boot (logit / inv.logit).
# =============================================================================

# Check and load required packages when this file is sourced.
for (pkg in c("mgcv", "VGAM", "boot")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' is required by simdatfun().")
}
library(mgcv); library(VGAM); library(boot)

# -----------------------------------------------------------------------------
# Arguments:
#   nclus     unused (kept in the signature so existing calls still work)
#   nspec     number of species (nspec == 10 reproduces the original exactly)
#   spnames   character vector of species short names, length nspec
#   nyr       number of years
#   p, disp   Tweedie power (p) and dispersion (disp = Tweedie phi / NB size)
#   qq        catchability matrix, ntactics rows x nspec columns
#   nvess     number of vessels
#   nlat,nlon spatial grid dimensions (need nlon > 1 for the E-W gradient models)
#   ntactics  number of fishing tactics (vessels cycle through 1:ntactics)
#   tripmode  central number of trips per vessel-year (actual = tripmode +/- 2)
#   spdist    spatial density model: "rndom", "ew_trend", or "ew_rndom"
#   rtype     observation distribution: "Tweedie" or "zinb"
# Returns a list: $sets (one row per set) and $B (nyr x nspec true biomass).
# Assumes total sets per vessel-year < 360 (holds for tripmode up to ~25).
# -----------------------------------------------------------------------------
simdatfun <- function(nclus, nspec, spnames, nyr, p, disp, qq,
                      nvess, nlat, nlon, ntactics, tripmode,
                      spdist = "rndom", rtype = "Tweedie") {
  year <- 1:nyr
  
  ## --- True biomass -------------------------------------------------------
  ## Each species follows a constant exponential trend: a small growth rate r
  ## of random sign, applied to a lognormal starting biomass B1. Abundant
  ## species start near 200, rarer ones near 50.
  r <- runif(nspec) * 0.1 * sample(c(-1, 1), nspec, replace = TRUE)  # rate in +/-[0,0.1]
  if (nspec == 10) {
    # Original ETBF mix: species 1-3 and 5 abundant (~200), the rest rarer (~50).
    B1 <- c(200 * exp(rnorm(3, 0, 0.5)), 50 * exp(rnorm(1, 0, 0.5)),
            200 * exp(rnorm(1, 0, 0.5)), 50 * exp(rnorm(5, 0, 0.5)))
  } else {
    # General rule: first ~40% of species abundant (~200), remainder rarer (~50).
    nbig <- max(1, round(0.4 * nspec))
    B1   <- ifelse(seq_len(nspec) <= nbig, 200, 50) * exp(rnorm(nspec, 0, 0.5))
  }
  B <- mat.or.vec(nyr, nspec)                         # nyr x nspec, filled below
  for (i in 1:nspec) B[, i] <- B1[i] * exp((year - 1) * r[i])
  
  ## --- Seasonal (monthly) catchability ------------------------------------
  ## qmon[month, species] scales catchability by month. The seasonal pattern is
  ## defined for the first 5 species; any further species get 1 (no season
  ## effect). Built as a 5 x 12 array, then transposed to 12 months x 5 species.
  qmon5 <- t(array(c(0.83, 0.49, 1.28, 0.98, 0.73,
                     0.82, 0.63, 1.17, 0.86, 0.46,
                     0.87, 0.92, 1.01, 1.38, 0.41,
                     0.95, 1.38, 0.72, 0.84, 0.58,
                     0.93, 1.81, 0.62, 0.71, 1.01,
                     0.97, 1.65, 0.63, 0.63, 1.71,
                     1.02, 1.45, 0.61, 0.59, 2.02,
                     1.56, 1.13, 0.85, 0.74, 1.72,
                     1.37, 0.79, 1.24, 0.81, 1.21,
                     0.99, 0.60, 1.32, 1.01, 0.79,
                     0.84, 0.50, 1.37, 1.55, 0.71,
                     0.85, 0.66, 1.19, 1.90, 0.65), dim = c(5, 12)))  # 12 mon x 5 sp
  if (nspec <= 5) {
    qmon <- as.data.frame(qmon5[, seq_len(nspec), drop = FALSE])
  } else {
    qmon <- as.data.frame(cbind(qmon5, matrix(1, nrow = 12, ncol = nspec - 5)))
  }
  colnames(qmon) <- spnames
  
  ## --- Vessels and their fishing tactics ----------------------------------
  ## Each vessel has a fixed base "tactic" (which species it targets, encoded by
  ## a row of qq). Vessels are assigned tactics 1,2,...,ntactics,1,2,... in turn.
  ## The lower half of the tactics "switch" to a higher tactic from year
  ## floor(nyr/2), mimicking a mid-series change in fleet behaviour. qq must have
  ## ntactics rows, so every tactic (original or switched) is a valid row index.
  shift     <- floor(ntactics / 2)                    # size of the tactic jump
  switchers <- seq_len(shift)                          # tactics that change (1..shift)
  vess <- data.frame(vess = 1:nvess)
  vess$tacst      <- rep(seq_len(ntactics), length.out = nvess)  # base tactic
  vess$tac_change <- NA                                # tactic after the switch
  vess$tcyr       <- NA                                # year the switch happens
  chg <- vess$tacst %in% switchers
  vess$tcyr[chg]       <- floor(nyr / 2)
  vess$tac_change[chg] <- vess$tacst[chg] + shift
  vess$vq <- exp(rnorm(nvess, mean = 0, sd = .3))      # per-vessel q: lognormal, median 1
  
  ## --- Vessel-years -------------------------------------------------------
  ## One row per (vessel, year). Draw the number of trips, copy the base tactic,
  ## then overwrite it with the switched tactic from year tcyr onwards.
  yrs    <- 1:nyr
  vessyr <- expand.grid(vess = vess$vess, yr = yrs)    # vess varies fastest
  vessyr$vykey <- paste(vessyr$vess, vessyr$yr, sep = "_")    # "vessel_year" key
  # Trips per vessel-year ~ uniform on tripmode +/- 2, floored at 1 so that
  # tripmode <= 2 cannot produce 0 (or negative) trips.
  vessyr$ntrip <- sample(seq(max(1, tripmode - 2), tripmode + 2),
                         size = nrow(vessyr), replace = TRUE)
  vessyr$tac   <- vess[vessyr$vess, ]$tacst            # base tactic per vessel-year
  av  <- vess[!is.na(vess$tac_change), ]$vess          # vessels that ever switch
  avy <- (1:nrow(vessyr))[vessyr$vess %in% av &        # rows at/after the switch year
                            vessyr$yr >= vess[vessyr$vess, ]$tcyr]
  vessyr$tac[avy] <- vess[vessyr$vess[avy], ]$tac_change
  
  ## --- Trips and sets -----------------------------------------------------
  ## Expand each vessel-year into its trips, and each trip into its sets. Built
  ## column-at-a-time (one data frame each) rather than row-by-row: rep()
  ## repeats each parent's fields, and sequence() numbers the children within
  ## each parent (e.g. sequence(c(3,2)) = 1,2,3,1,2). Keys (vykey, tripkey) let
  ## us join the levels back together later.
  ntr  <- vessyr$ntrip                                 # trips per vessel-year
  trip <- data.frame(vess  = rep(vessyr$vess,  ntr),
                     yr    = rep(vessyr$yr,    ntr),
                     trip  = sequence(ntr),            # 1..ntrip within each vessel-year
                     vykey = rep(vessyr$vykey, ntr),
                     nsets = sample(5:13, sum(ntr), replace = TRUE))
  trip$tripkey <- paste(trip$vess, trip$yr, trip$trip, sep = "_")   # "vessel_year_trip"
  
  ns   <- trip$nsets                                   # sets per trip
  sets <- data.frame(vess    = rep(trip$vess,    ns),
                     yr      = rep(trip$yr,      ns),
                     trip    = rep(trip$trip,    ns),
                     tripkey = rep(trip$tripkey, ns),
                     set     = sequence(ns))           # 1..nsets within each trip
  
  ## --- Place trips/sets on a calendar -------------------------------------
  ## Give each set a day-of-year (1..360). Total sets per vessel-year are summed
  ## (joined back by key, so row order can't matter), leaving 360 - nsets "spare"
  ## days. Trip start days are sampled within that window; a trip's sets then run
  ## on consecutive days.
  a <- aggregate(nsets ~ vykey, data = trip, FUN = sum)
  vessyr$nsets   <- a$nsets[match(vessyr$vykey, a$vykey)]
  vessyr$shortyr <- 360 - vessyr$nsets
  
  # Sample ntrip distinct, sorted start-day offsets per vessel-year.
  # lapply (not sapply) guarantees a list even if all ntrip happen to be equal.
  mkdates   <- function(x) sort(sample(1:vessyr$shortyr[x], size = vessyr$ntrip[x],
                                       replace = FALSE))
  tripdates <- lapply(1:nrow(vessyr), mkdates)
  
  # Turn each vessel-year's start offsets + trip lengths into the first (tripst)
  # and last (tripnd) day of every trip; b1 staggers trips so they don't overlap.
  tripstart <- function(x) {
    a  <- tripdates[[x]]
    b  <- trip[trip$vykey == vessyr$vykey[x], "nsets"]   # sets per trip, in trip order
    b1 <- cumsum(b - 1)
    cbind(tripst = 1 + a + b1 - b, tripnd = a + b1)
  }
  trip <- cbind(trip, do.call(rbind, lapply(1:nrow(vessyr), tripstart)))
  
  # Expand each trip's [tripst, tripnd] into one day per set, in trip order.
  # Month = day/30, clamped to 12 so it always indexes the 12-row qmon table.
  dayofset  <- function(x) trip[x, "tripst"]:trip[x, "tripnd"]
  sets$day  <- do.call(c, lapply(1:nrow(trip), dayofset))
  sets$mon  <- pmin(floor(sets$day / 30) + 1, 12)
  
  ## --- Tactic per set -----------------------------------------------------
  ## Join the vessel-year tactic onto trips (by vykey), then onto sets (by
  ## tripkey). match() looks the key up BY VALUE - indexing a vector with the
  ## key string directly would misalign (the original bug).
  trip$tactic <- vessyr$tac[match(trip$vykey, vessyr$vykey)]
  sets$tactic <- trip$tactic[match(sets$tripkey, trip$tripkey)]
  
  ## --- Spatial density per grid cell --------------------------------------
  ## relden has one row per (lat, lon) cell and one column per species, giving
  ## each species' relative density there. expand.grid varies lat fastest.
  relden <- expand.grid(lat = 1:nlat, lon = 1:nlon)
  relden$loc <- with(relden, paste(lat, lon, sep = "_"))    # "lat_lon" cell key
  
  # "rndom": density is uniform noise (no spatial structure).
  do_rndom <- function(rmin = 0, rmax = 1)
    matrix(runif(nspec * nlat * nlon, min = rmin, max = rmax),
           ncol = nspec, nrow = nrow(relden))
  
  # "ew_trend": each species has a linear west->east gradient between random
  # west and east endpoints. Density varies with LONGITUDE, constant across
  # latitude. Row order matches relden; correct for any nlat, nlon (need nlon>1).
  do_ewtrend <- function() {
    west <- runif(nspec); east <- runif(nspec)
    grad <- (east - west) / (nlon - 1)                 # E-W increment per lon step
    outer(relden$lon - 1, grad) +                      # cell value = (lon-1)*grad ...
      matrix(west, nrow = nrow(relden), ncol = nspec, byrow = TRUE)  # ... + west
  }
  
  # Choose the spatial model. "ew_rndom" blends a gradient with noise on the
  # logit scale, which keeps the result in (0,1).
  if (spdist == "rndom")    lldist <- do_rndom()
  if (spdist == "ew_trend") lldist <- do_ewtrend()
  if (spdist == "ew_rndom") lldist <- inv.logit(logit(do_ewtrend()) + .2 * logit(do_rndom(0, 1)))
  
  # Normalise each species' density to mean 1 across cells. This sets the overall
  # catch scale (and makes it grid-size independent) WITHOUT altering the pattern.
  lldist <- sweep(lldist, 2, colMeans(lldist), "/")
  relden <- cbind(relden, lldist)
  spcols <- 3 + seq_len(nspec)                         # density columns (after lat,lon,loc)
  
  ## --- Place each set in space --------------------------------------------
  ## The first set of every trip starts at a random cell; later sets follow a
  ## small random walk from the previous set, so a trip stays roughly local.
  ## Done on plain integer vectors (not data-frame rows) for speed: all step
  ## increments are drawn up front, then a cheap sequential pass applies them
  ## with reflection at the grid edges. The pass stays a loop because the clamp
  ## needs the running position (so it can't be a plain cumsum).
  n     <- nrow(sets)
  first <- sets$set == 1                               # first set of each trip
  lat <- integer(n); lon <- integer(n)
  lat[first] <- sample(1:nlat, sum(first), replace = TRUE)
  lon[first] <- sample(1:nlon, sum(first), replace = TRUE)
  dlat <- sample(c(-1L, 0L, 1L), n, replace = TRUE, prob = c(.05, .9, .05))
  dlon <- sample(c(-1L, 0L, 1L), n, replace = TRUE, prob = c(.05, .9, .05))
  for (j in which(!first)) {                           # j-1 = previous set, same trip
    np <- lat[j - 1] + dlat[j]; lat[j] <- if (np < 1 || np > nlat) lat[j - 1] else np
    np <- lon[j - 1] + dlon[j]; lon[j] <- if (np < 1 || np > nlon) lon[j - 1] else np
  }
  sets$lat <- lat; sets$lon <- lon
  sets$loc <- paste(lat, lon, sep = "_")
  
  ## --- Catch --------------------------------------------------------------
  ## Realised catchability per set = tactic q x vessel q x (1) x month q. The
  ## middle "1" is a placeholder for any location-specific q. The vessel-q
  ## vector recycles down the columns, scaling each row by that set's vessel.
  qsp <- qq[sets$tactic, ] * vess$vq[sets$vess] * 1 * qmon[sets$mon, ]
  
  # Expected catch per set/species = effort(=1) x biomass(year) x density(cell) x q.
  catch_M <- 1 * B[sets$yr, ] * relden[match(sets$loc, relden$loc), spcols] * qsp
  colnames(catch_M) <- paste0(spnames, "_M")           # "_M" = expected (Mean)
  
  # Observed catch: draw from the chosen distribution, one species column at a time.
  gendat <- function(x) switch(rtype,
                               "Tweedie" = rTweedie(x, p = p, phi = disp),
                               "zinb"    = rzinegbin(n = length(x), size = disp, mu = x, pstr0 = 0.05))
  catchobs <- round(apply(catch_M, 2, gendat), 0)
  colnames(catchobs) <- paste0(spnames, "_R")          # "_R" = realised/observed
  
  # Keep the set-description columns (through "loc") and append the catches.
  a <- grep("loc", names(sets))
  sets <- cbind(sets[, 1:a], catch_M, catchobs)
  return(list(sets = sets, B = B))
}