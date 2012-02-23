################################################################################
## Wrap shogun static with some passible R interface

## Some info
## Two class labels should be real/numeric -1,1
## multiclass labels should be real/numeric [0, 1, 2, ..]
## features of (obs x features) are sent into `sg` as their transpose!
.set.loglevel <- function(verbose) {
  level <- if (verbose) 'INFO' else 'INFO'
  sg('loglevel', level)
}

matchShogunKernel <- function(kernel) {
  kernels <- c('spectrum', 'weighted-degree', 'WD',
               'weighted-degree-with.shifs', 'WDS',
               'gaussian', 'polynomial', 'sigmoid')
}

setGeneric("SVM", signature=c("x", "kernel"),
function(x, kernel, C=1, C2=C, nu=0.2, epsilon=0.1, class.weights=NULL,
         cross=0, cache=40, tol=0.001, shrinking=TRUE, ..., verbose=FALSE) {
  standardGeneric("SVM")
})

setMethod("SVM", c(x="DNAStringSet"),
function(x, kernel, C..., verbose=FALSE) {
  .set.loglevel(verbose)

})

setMethod("SVM", c(x="matrix"),
function(x, kernel, ..., verbose) {

})

setMethod("predict", c(object="SVM"),
function(object, newdata, ...) {

})

SVM <- function(x, kernel=...)

################################################################################
## Hide me

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
