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


    # read in dates
    dat <- read.csv(paste0(top.dir.path, i, "/", i, ".csv"))
    
    # read in parameters
    pars <- read.csv(paste0(top.dir.path, i, "/", "bacon_pars", ".csv"))

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



#' Title
#'
#' @param top.dir.path 
#'
#' @return
#' @export
#'
#' @examples
#' summary.age.mods <- AggregateSummaryAgeModels("inst/extdata/terr_14C_min10_dates-2020.03.04_15-19-42/") %>%
#'   tbl_df()
#' 
#' summary.age.mods %>%
#'   filter(DataName %in% sample(unique(DataName), 5)) %>%
#'   ggplot(aes(x = depth, y = mean, colour = DataName)) +
#'   geom_ribbon(aes(ymax = max, ymin = min, fill = DataName), colour = NA, alpha = 0.25) +
#'   geom_line() +
#'   geom_line(aes(y = median), linetype = 2)
AggregateSummaryAgeModels <- function(top.dir.path){
  
  # Get list of directories inside top-level dir
  dirs <- list.dirs(top.dir.path, full.names = FALSE, recursive = FALSE)
  paths <- list.dirs(top.dir.path, full.names = TRUE, recursive = FALSE)
  
  #browser()
  
  fls <- lapply(paths, function(i) {
    
    ages.flnm <- Sys.glob(file.path(i, "*ages.txt"))
    #print(ages.flnm)

    if (identical(ages.flnm, character(0)) == FALSE){
      ages.df <- read.delim(ages.flnm)
      
      pars.df <- read.csv(file.path(i, "bacon_pars.csv"))
      
      cbind(pars.df, ages.df)
    }else if (identical(ages.flnm, character(0))){
      
      df.null <- data.frame(
        depth = NA, min = NA, max = NA, median = NA, mean = NA
      )
      pars.df <- read.csv(file.path(i, "bacon_pars.csv"))
      
      cbind(pars.df, df.null)
      
      }
    
  })
  df <- do.call(rbind.data.frame, fls)
  
  return(df)
}


.StackIterations <- function(x){
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  n.row <- nrow(x)
 
  x <- stack(x, stringsAsFactors = FALSE)
  x$depth.index <- (1:n.row) -1
  x
}

#' Title
#'
#' @param bacon.posterior 
#' @param pars 
#'
#' @return
#' @export
#'
#' @examples
#' posterior <- read.table("inst/extdata/test-subset/AE3/AE3_62.out", header = FALSE)
#' pars <- read.csv("inst/extdata/test-subset/AE3/bacon_pars.csv", stringsAsFactors = FALSE)
#' 
#' 
#' age.mods <- ConstructAgeModels(posterior, pars) %>%
#'   as_tibble()
#' 
#' age.mods %>%
#'   filter(iter %in% sample(unique(iter), 100)) %>%
#'   ggplot(aes(x = depth, y = age, group = iter)) +
#'   geom_line(alpha = 0.1) +
#'   expand_limits(x = 0, y = 0)
ConstructAgeModels <- function(bacon.posterior, pars){
  
  n.col <- ncol(bacon.posterior)
  log.lik <- bacon.posterior[, n.col]
  
  
  bacon.posterior[, 2:(n.col - 2)] <- bacon.posterior[, 2:(n.col - 2)] * pars$thick
  
  age.mods <- apply(bacon.posterior[, 1:(n.col - 2)], 1, cumsum)
 
  age.mods <- .StackIterations(age.mods)
  
  age.mods$depth <- pars$d.min + age.mods$depth.index * pars$thick
  
  age.mods <- age.mods[, c("ind", "depth", "values")]
  names(age.mods) <- c("iter", "depth", "age") 
  
  return(age.mods)
  
  }



AggregatePosteriorAgeModels <- function(top.dir.path){
  
  # Get list of directories inside top-level dir
  dirs <- list.dirs(top.dir.path, full.names = FALSE, recursive = FALSE)
  paths <- list.dirs(top.dir.path, full.names = TRUE, recursive = FALSE)
  
  #browser()
  
  fls <- lapply(paths, function(i) {
    
    posterior.flnm <- Sys.glob(file.path(i, "*.out"))
    #print(ages.flnm)
    
    if (identical(posterior.flnm, character(0)) == FALSE){
      posterior <- read.table(posterior.flnm, header = FALSE)
      
      pars.df <- read.csv(file.path(i, "bacon_pars.csv"))
      
      age.mods <- ConstructAgeModels(posterior, pars.df)
      cbind(pars, age.mods)
      
    }else if (identical(posterior.flnm, character(0))){
      
      df.null <- data.frame(
        iter = NA, depth = NA, age = NA
      )
      pars.df <- read.csv(file.path(i, "bacon_pars.csv"))
      
      cbind(pars.df, df.null)
      
    }
    
  })
  df <- do.call(rbind.data.frame, fls)
  
  return(df)
}


