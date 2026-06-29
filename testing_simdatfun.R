source("simdatfun.R", keep.source = TRUE)   # keep.source enables line profiling
set.seed(42)
spnames <- c("ALB","BET","YFT","SWO","MLS","DOL","SBT","BUM","BLM","SHK")
qq <- rbind(c(0.50,1.00,0.80,0.05,0.01,0.01,0.01,0.01,0.01,0.30),
            c(0.20,0.10,0.50,1.00,0.50,0.01,0.01,0.60,0.60,0.60),
            c(1.00,0.40,1.00,0.10,1.00,0.05,0.20,0.20,0.30,0.30),
            c(0.70,0.50,0.01,0.10,0.10,0.50,1.00,0.20,0.20,0.60))

run1 <- function()
  simdatfun(4, 10, spnames, nyr=20, p=1.3, disp=10, qq=qq, nvess=16,
            nlat=10, nlon=10, ntactics=4, tripmode=10,
            spdist="ew_rndom", rtype="Tweedie")

library(profvis)
p <- profvis({ for (i in 1:10) run1() }, interval = 0.005)
p                                   # interactive flame graph + per-line table
# htmlwidgets::saveWidget(p, "profile.html")   # save to view/share

p <- profvis({ for (i in 1:10) run1() })
writeLines(p$x$message$prof_output)              # may be NULL depending on version
# More reliable — re-run with Rprof to get text directly:
tmp <- tempfile()
Rprof(tmp, interval = 0.005, line.profiling = TRUE)
for (i in 1:10) invisible(run1())
Rprof(NULL)
s <- summaryRprof(tmp, lines = "show")
print(head(s$by.self, 15))
print(head(s$by.line, 20))