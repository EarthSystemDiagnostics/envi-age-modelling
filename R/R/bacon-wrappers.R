## Batch functions for Bacon ---------

#' Create folder and required files for each site
#' @description Creates a set of directories and the files required to run
#' rbacon::Bacon taking data for multiple sites from a dataframe.
#' @param dat A dataframe
#' @param filename
#' @param path
#' @param suffix Append a suffix to the top-level folder name. Defaults to the
#'  current date. Alternatives are current date and time "date.time" or "none".
#' @param site.id Name of variable uniquily identifying sites
#' @param sample.id Optional sample.id for Bacon
#' @param age Name of variable containing the age estimates
#' @param age.err Name of variable containing the 1-sigma age uncertainties
#' @param depth Name of variable containing the depth of the samples
#'
#' @return NULL - has side effects of creating a folder structure.
#' @export
#'
#' @examples
MakeBaconDirs <- function(dat = NULL, filename, path,
                              suffix = c("date", "date.time", "none"),
                              site.id, sample.id, age, age.err, depth){

  if (is.null(dat)){
    dat <- read.table(paste0(path, filename),
                      header = TRUE,
                      sep = ";",
                      dec = ",",
                      stringsAsFactors = FALSE
    )
  }else{
    dat <- dat
  }


  # split data.frame into a list, each element is one site
  dat.lst <- split(dat, dat[[site.id]])

  # names of sites
  nms <- names(dat.lst)

  # create bacon compatible dataframe for each
  bacon.dat.lst <- lapply(dat.lst, function(i) {

    d <- data.frame(
      id = i[[sample.id]],
      age = i[[age]],
      error = i[[age.err]],
      depth = i[[depth]]
    )

    d <- subset(d, complete.cases(d))

    d
  })


  # create top-level folder name

  suffix <- match.arg(suffix)
  suff <- switch(suffix,
                 none = "",
                 date = paste0("-", Sys.Date()),
                 date.time = paste0("-", format(Sys.time(), "%Y.%m.%d_%H-%M-%S"))
  )

  folder.name <- paste0(
    # remove suffice from filename
    strsplit(filename, split = ".", fixed = TRUE)[[1]][1],
    # add current date and time
    #"-", format(Sys.time(), "%Y.%m.%d_%H-%M-%S")
    suff
  )

  # create top-level directory named data.file name
  folder.path <- paste0(path, "/", folder.name)
  dir.create(folder.path)


  folder.paths <- lapply(nms, function(i) paste0(folder.path, "/", i))

  lapply(folder.paths, function(i) dir.create(path = i))

  # write bacon input files
  lapply(1:length(bacon.dat.lst), function(i) {
    write.csv(bacon.dat.lst[[i]],
              file = paste0(folder.paths[i], "/", nms[i], ".csv"),
              row.names = FALSE, quote = FALSE)
  })
}


CreateParametersFiles <- function(top.dir.path, pars.df){

  # Get list of directories inside top-level dir
  dirs <- list.dirs(top.dir.path, full.names = FALSE, recursive = FALSE)

  dir.paths <- list.dirs(top.dir.path, full.names = TRUE, recursive = FALSE)

  pars.present <- dirs %in% pars.df$DataName

  if (any(pars.present == FALSE))  stop(paste0("parameters missing for ", dirs[pars.present == FALSE]))


  lapply(1:length(dirs), function(i){

    pars <- subset(pars.df, pars.df$DataName == dirs[[i]])

    write.csv(pars,
              file = paste0(dir.paths[i], "/", "bacon_pars", ".csv"),
              row.names = FALSE, quote = FALSE)

    })

}




#' Run rbacon::Bacon for data in a set of directories
#'
#' @param top.dir.path Top level directory under which are multiple directories
#' containing input files for Bacon
#'
#' @return NULL - has side effects of writing Bacon output to a folder structure.
#' @export
#'
#' @examples
RunBaconDirs <- function(top.dir.path, runname = "", frac.cores = 0.5){

  if (.Platform$OS.type == "unix") {
    n.cores.available <- parallel::detectCores()
    n.cores <- ceiling(n.cores.available * frac.cores)
    }else{n.cores <- 1}

  # Get list of directories inside top-level dir
  dirs <- list.dirs(top.dir.path, full.names = FALSE, recursive = FALSE)

  out <- parallel::mclapply(dirs, function(i) {


    # read in the data so that Bacon parameters can be calculated

    # we can set more parameters here, like the depths to be interpolated at the
    # end and the start and end depths

    dat <- read.csv(paste0(top.dir.path, i, "/", i, ".csv"))

    pars <- read.csv(paste0(top.dir.path, i, "/", "bacon_pars", ".csv"))

    # # estimate mean sediment accumulation rate
    # # use simple linear regression
    # lm1 <- lm(age~depth, data = dat)
    # acc.mean <- coef(lm1)[2]
    #
    # # calculate mean depth step and total core length to use for
    # # setting a good layer thickness
    # geo.mean.d.depth <- (mean(sqrt(diff(sort(dat$depth)))))^2
    # core.length <- diff(range(dat$depth))
    #
    # # set layer thickness to a fraction of median depth step
    # thick = geo.mean.d.depth / 3

    # call Bacon
    Bacon2(core = i, coredir = top.dir.path,
           d.min = pars$d.min, d.max = pars$d.max, d.by = pars$d.by,
           acc.mean = pars$acc.mean, thick = pars$thick,
           # suppress interactive questions
           suggest = FALSE, ask = FALSE,
           runname = runname, remember = FALSE,
           plot.pdf = TRUE, suppress.plots = TRUE,
           verbose = FALSE)

    return(i)

  }, mc.cores = n.cores, mc.preschedule = FALSE)
  return(out)
}
