#' Bacon2, A Wrapper for rbacon::Bacon so that it plays nicely when run in
#' parallel. Suppresses calls to graphics devices.
#'
#' @param suppress.plots If TRUE it prevents calls to on screen graphics devices.
#' @inheritParams rbacon::Bacon
#' @param ...
#'
#' @return
#' @import rbacon
#' @export
#'
#' @examples
Bacon2 <- function (suppress.plots = TRUE,
                    core = "MSB2K", thick = 5, coredir = "", prob = 0.95,
                    d.min = NA, d.max = NA, d.by = 1, seed = NA, depths.file = FALSE,
                    depths = c(), depth.unit = "cm", age.unit = "yr", unit = depth.unit,
                    acc.shape = 1.5, acc.mean = 20, mem.strength = 4, mem.mean = 0.7,
                    boundary = NA, hiatus.depths = NA, hiatus.max = 10000, add = c(),
                    after = 1e-04/thick, cc = 1, cc1 = "IntCal13", cc2 = "Marine13",
                    cc3 = "SHCal13", cc4 = "ConstCal", ccdir = "", postbomb = 0,
                    delta.R = 0, delta.STD = 0, t.a = 3, t.b = 4, normal = FALSE,
                    suggest = TRUE, reswarn = c(10, 200), remember = TRUE, ask = TRUE,
                    run = TRUE, defaults = "defaultBacon_settings.txt", sep = ",",
                    dec = ".", runname = "", slump = c(), BCAD = FALSE, ssize = 2000,
                    th0 = c(), burnin = min(500, ssize), MinAge = c(), MaxAge = c(),
                    MinYr = MinAge, MaxYr = MaxAge, cutoff = 0.001, plot.pdf = TRUE,
                    dark = 1, date.res = 100, age.res = 200, yr.res = age.res,
                    close.connections = TRUE, verbose = TRUE, ...)
{

  # Temporarily attach all internal functions from package rbacon
  attach(loadNamespace("rbacon"), name = "rbacon_all")

  coredir <- assign_coredir(coredir, core, ask)
  if (core == "MSB2K" || core == "RLGH3") {
    dir.create(paste(coredir, core, "/", sep = ""), showWarnings = FALSE,
               recursive = TRUE)
    fileCopy <- system.file(paste("extdata/Cores/", core,
                                  sep = ""), package = "rbacon")
    file.copy(fileCopy, coredir, recursive = TRUE, overwrite = FALSE)
  }
  if (ccdir == "")
    ccdir <- paste(system.file("extdata", package = "rbacon"),
                   "/Curves/", sep = "")
  ccdir <- .validateDirectoryName(ccdir)
  defaults <- system.file("extdata", defaults, package = "rbacon")
  dets <- .read.dets(core, coredir, sep = sep, dec = dec,
                     cc = cc)
  if (ncol(dets) > 4 && length(cc) > 0) {
    cc.csv <- unique(dets[, 5])
    if (verbose) {
      if (length(cc.csv) == 1) {
        if (cc.csv != cc)
          message(" Using calibration curve specified within the .csv file,",
                  cc[cc.csv], "\n")
      }
      else if (min(cc.csv) == 0)
        message(" Using a mix of cal BP and calibrated C-14 dates\n")
      else message(" Using several C-14 calibration curves\n")
    }
  }
  if (suggest) {
    sugg <- sapply(c(1, 2, 5), function(x) x * 10^(-1:2))
    ballpacc <- lm(dets[, 2] * 1.1 ~ dets[, 4])$coefficients[2]
    ballpacc <- abs(sugg - ballpacc)
    ballpacc <- ballpacc[ballpacc > 0]
    sugg <- sugg[order(ballpacc)[1]]
    if (!sugg %in% acc.mean) {
      ans <- readline(message(" Ballpark estimates suggest changing the prior for acc.mean to ",
                              sugg, " ", age.unit, "/", depth.unit, ". OK? (y/N) "))
      if (tolower(substr(ans, 1, 1)) == "y")
        acc.mean <- sugg
      else message(" No problem, using the provided prior")
    }
  }
  if (!is.na(boundary[1]))
    boundary <- sort(unique(boundary))
  if (!is.na(hiatus.depths[1])) {
    hiatus.depths <- sort(unique(hiatus.depths))
    if (length(acc.mean) == 1)
      acc.mean <- rep(acc.mean, length(hiatus.depths) +
                        1)
  }
  info <- .Bacon.settings(core = core, coredir = coredir,
                          dets = dets, thick = thick, remember = remember, d.min = d.min,
                          d.max = d.max, d.by = d.by, depths.file = depths.file,
                          slump = slump, acc.mean = acc.mean, acc.shape = acc.shape,
                          mem.mean = mem.mean, mem.strength = mem.strength, boundary = boundary,
                          hiatus.depths = hiatus.depths, hiatus.max = hiatus.max,
                          BCAD = BCAD, cc = cc, postbomb = postbomb, cc1 = cc1,
                          cc2 = cc2, cc3 = cc3, cc4 = cc4, depth.unit = depth.unit,
                          normal = normal, t.a = t.a, t.b = t.b, delta.R = delta.R,
                          delta.STD = delta.STD, prob = prob, defaults = defaults,
                          runname = runname, ssize = ssize, dark = dark, MinAge = MinAge,
                          MaxAge = MaxAge, cutoff = cutoff, age.res = age.res,
                          after = after, age.unit = age.unit)
  .assign_to_global("info", info)
  info$coredir <- coredir
  info$seed <- seed
  info$isplum <- FALSE
  if (any(info$acc.shape == info$acc.mean))
    stop("acc.shape cannot be equal to acc.mean", call. = FALSE)
  if (info$t.b - info$t.a != 1)
    stop("t.b - t.a should always be 1, check the manual",
         call. = FALSE)
  if (min(acc.shape) < 1)
    cat("\nWarning, using values <1 for acc.shape might cause unexpected results\n")
  if (info$cc > 0)
    if (info$postbomb == 0 && ((ncol(info$dets) == 4 &&
                                min(info$dets[, 2]) < 0) || ncol(info$dets) > 4 &&
                               max(info$dets[, 5]) > 0 && min(info$dets[info$dets[,
                                                                                  5] > 0, 2]) < 0))
      stop("you have negative C14 ages so should select a postbomb curve",
           call. = FALSE)
  info$calib <- .bacon.calib(dets, info, date.res, ccdir = ccdir)
  info$rng <- c()
  for (i in 1:length(info$calib$probs)) {
    tmp <- info$calib$probs[[i]]
    info$rng <- range(info$rng, tmp[which(tmp[, 2] > cutoff),
                                    1])
  }
  if (length(th0) == 0)
    info$th0 <- round(rnorm(2, max(MinAge, dets[1, 2]),
                            dets[1, 3]))
  info$th0[info$th0 < info$MinAge] <- info$MinAge
  if (length(depths) == 0)
    depths <- seq(info$d.min, info$d.max, by = d.by)
  if (depths.file) {
    dfile <- paste0(info$coredir, info$core, "/", info$core,
                    "_depths.txt")
    if (!file.exists(dfile))
      stop("I cannot find the file ", paste0(info$coredir,
                                             info$core, "/", info$core, "_depths.txt"), call. = FALSE)
    depths <- read.table(dfile, header = FALSE)[, 1]
    if (!is.numeric(depths[1]))
      stop("File should contain numbers only, no headers",
           call. = FALSE)
  }
  info$depths <- depths
  if (min(depths) < info$d.min)
    info$d.min <- min(depths)
  if (max(depths) > info$d.max)
    info$d.max <- max(depths)
  info$elbows <- seq(floor(info$d.min), ceiling(info$d.max),
                     by = thick)
  info$K <- length(info$elbows)
  info$cK <- info$d.min + (info$thick * info$K)
  if (!is.na(info$hiatus.depths[1]) || !is.na(info$boundary[1])) {
    ifelse(is.na(info$boundary[1]), hd <- info$hiatus.depths,
           hd <- info$boundary)
    if (min(hd) < info$d.min)
      stop("cannot have hiatus above the core's top depth. Adapt hiatus.depths or d.min.",
           call. = FALSE)
    if (max(hd) + info$thick > info$d.max)
      stop("the age-depth model should have at least one section below the one containing the deepest hiatus. Adapt thick or d.max?",
           call. = FALSE)
    if (length(hd) > 1) {
      above <- c()
      for (i in hd) above <- c(above, max(which(info$elbows <=
                                                  i)))
      if (any(diff(above) < 2))
        stop("we need at least 2 section elbows between hiatuses. Choose fewer hiatuses, different depths, more sections (decrease thick) or a different d.min.\n ",
             call. = FALSE)
    }
  }
  ans <- "n"
  if (suggest)
    if (length(reswarn) == 2)
      if (info$K < min(reswarn)) {
        sugg <- pretty(thick * (info$K/min(reswarn)),
                       10)
        sugg <- min(sugg[sugg > 0])
        ans <- readline(message(" Warning, the current value for thick, ",
                                thick, ", will result in very few age-model sections (",
                                info$K, ", not very flexible). Suggested maximum value for thick: ",
                                sugg, " OK? (y/n) "))
      }
  else if (info$K > max(reswarn)) {
    sugg <- max(pretty(thick * (info$K/max(reswarn))))
    ans <- readline(message(" Warning, the current value for thick, ",
                            thick, ", will result in very many age-model sections (",
                            info$K, ", possibly hard to run). Suggested minimum value for thick: ",
                            sugg, " OK? (y/n) "))
  }
  if (tolower(substr(ans, 1, 1)) == "y") {
    cat(" OK, setting thick to ", sugg, "\n")
    thick <- sugg
    info$thick = thick
    info$elbows <- seq(floor(info$d.min), ceiling(info$d.max),
                       by = thick)
    if (length(info$slump) > 0)
      info$elbows <- seq(floor(info$d.min), toslump(ceiling(info$d.max),
                                                    info$slump), by = thick)
    info$K <- length(info$elbows)
    info$cK <- info$d.min + (info$thick * info$K)
  }
  if (length(slump) > 0) {
    if (length(slump)%%2 == 1)
      stop("slumps need both upper and lower depths. Please check the manual",
           call. = FALSE)
    slump <- matrix(sort(slump), ncol = 2, byrow = TRUE)
    info$slump <- slump
    slumpdmax <- toslump(ceiling(info$d.max), slump)
    info$elbows <- seq(floor(info$d.min), slumpdmax, by = thick)
    info$K <- length(info$elbows)
    info$cK <- info$d.min + (info$thick * info$K)
    info$slumpfree <- toslump(depths, slump)
    info$slumphiatus <- toslump(info$hiatus.depths, slump)
    if (!is.na(info$boundary[1])) {
      info$slumpboundary <- toslump(info$boundary, slump)
      info$slumphiatus <- info$slumpboundary
    }
    slumpdets <- info$dets
    slumpdets[, 4] <- toslump(slumpdets[, 4], slump, remove = FALSE)
    info$slumpdets <- slumpdets[!is.na(slumpdets[, 4]),
                                ]
  }
  info$prefix <- paste(coredir, core, "/", core, runname,
                       "_", info$K, sep = "")
  info$coredir <- coredir
  info$bacon.file <- paste(info$prefix, ".bacon", sep = "")
  if (!file.exists(outfile <- paste(info$prefix, ".out", sep = "")))
    file.create(outfile)
  if (BCAD)
    info$BCAD <- TRUE
  if (!is.na(boundary[1])) {
    if (length(slump) > 0)
      boundary <- info$slumpboundary
    info$hiatus.depths <- boundary
    if (length(add) == 0)
      add <- info$acc.mean
    info$hiatus.max <- add
  }
  .assign_to_global("info", info)
  prepare <- function() {
    if (suppress.plots == FALSE){
      pn <- c(1, 2, 3, 3)
      if (!is.na(info$hiatus.depths[1]))
        if (is.na(info$boundary[1]))
          pn <- c(1, 2, 3, 4, 4, 4)
        layout(matrix(pn, nrow = 2, byrow = TRUE), heights = c(0.3,
                                                               0.7))
        oldpar <- par(mar = c(3, 3, 1, 1), mgp = c(1.5, 0.7,
                                                   0), bty = "l")
        on.exit(par(oldpar))
        .PlotAccPrior(info$acc.shape, info$acc.mean, depth.unit = depth.unit,
                      age.unit = age.unit)
        .PlotMemPrior(info$mem.strength, info$mem.mean, thick)
        if (!is.na(info$hiatus.depths)[1])
          if (is.na(info$boundary)[1])
            .PlotHiatusPrior(info$hiatus.max, info$hiatus.depths)
        calib.plot(info, BCAD = BCAD)
        legend("topleft", core, bty = "n", cex = 1.5)
    }
  }
  cook <- function() {
    txt <- paste(info$prefix, ".bacon", sep = "")
    bacon(txt, outfile, ssize, ccdir)
    scissors(burnin, info)
    if (suppress.plots == FALSE){
      agedepth(info, BCAD = BCAD, depths.file = depths.file,
               depths = depths, verbose = TRUE, age.unit = age.unit,
               depth.unit = depth.unit, ...)
    }
    if (plot.pdf){
      pdf(file = paste0(info$prefix, ".pdf"))
      agedepth(info, BCAD = BCAD, depths.file = depths.file,
               depths = depths, verbose = FALSE, age.unit = age.unit,
               depth.unit = depth.unit, ...)
      dev.off()
    }
  }
  .write.Bacon.file(info)
  if (!run)
    prepare()
  else if (!ask)
    cook()
  else {
    prepare()
    ans <- readline(message(" Run ", core, " with ", info$K,
                            " sections? (Y/n) "))
    ans <- tolower(substr(ans, 1, 1))[1]
    if (ans == "y" || ans == "")
      cook()
    else message("  OK. Please adapt settings")
  }
  if (close.connections)
    closeAllConnections()

  # detach internal rbacon runctions
  detach("rbacon_all")
}
