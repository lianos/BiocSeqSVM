##' Fix shogun r_static installation.
##' 
##' The installation of the shogun static library doesn't work
##' with new versions of R.
##' 
##' First: the \code{sg/Meta/package.rds} file is broken. A 'Built' entry
##' needs to be added to its $DESCRIPTION vector.
##' 
##' Second: The compiled ligrary \code{sg.so} isn't save in its architture-
##' specific subfolder. For example, on OSX, the \code{sg.so} file is saved in
##' \code{sg/libs}, and needs to be moved to \code{sg/libs/{x86_64|i386}}.
fix.shogun.static.install <- function() {
  cat("... fixing package.rds\n")
  ## The `package.rds` file as compiled is hosed, need to add
  ## a 'built' string to the `$DESCRIPTION` character vector, otherwise
  ## installing new packages is b0rked, eg:
  ## 
  ##   R> install.packages('Matrix')
  ##   Error in .readPkgDesc(lib, fields) : 
  ##     number of items to replace is not a multiple of replacement length
  rds.path <- system.file('Meta', 'package.rds', package="sg")
  if (!file.exists(rds.path)) {
    stop("package.rds file not found:\n  ", rds.path)
  }
  x <- readRDS(rds.path)
  x$DESCRIPTION['Built'] <- paste(R.version$major, R.version$minor, sep=".")
  saveRDS(x, rds.path)
  
  
  cat("... moving sg.so to appropriate place\n")
  r.arch <- .Platform$r_arch
  bad.so.path <- system.file('libs', 'sg.so', package="sg")
  if (!file.exists(bad.so.path)) {
    cat("sg.so not found, did you already fix this?\n")
  } else {
    lib.dir <- dirname(bad.so.path)
    cmd <- paste('file', bad.so.path)
    info <- system(cmd, intern=TRUE)
    arch.match <- regexpr(r.arch, info)
    if (arch.match > -1) {
      arch.dir <- file.path(lib.dir, r.arch)
      if (!dir.exists(arch.dir)) {
        if (dir.create(arch.dir)) {
          
        } else {
          warning("Could not create library directory ", arch.dir)
        }
      }
    } else {
      warning("sg.so architecture does not seem to match R_arch: ", r.arch,
              "\n")
    }
  }
  
}
