# This is the R part of the irlba package embedded into one file, based on its
# state in Github at https://github.com/bwlewis/irlba/tree/bb71306fc3c17e4325b8e80aa627ee60c8edde01
# It is licensed under the GPL-3. See https://cran.r-project.org/package=irlba
# for more details and the authors and contributors.
# Due to persistent issues where a change to the Matrix package broke irlba
# (see https://github.com/bwlewis/irlba/issues/70), I am taking the nuclear
# options and embedding the package into the code base. This should fix the
# problem, at the cost of losing the C-based fast code path.
# Apart from removing the fastpath code, I have made only minimal cosmetic
# changes, including renaming functions to make them easier to search for if I
# am ever able to reverse this change, and also no longer exporting any of the
# functions.

# Find a few approximate singular values and corresponding
# singular vectors of a matrix.
#
# The augmented implicitly restarted Lanczos bidiagonalization algorithm
# (IRLBA) finds a few approximate largest (or, optionally, smallest) singular
# values and corresponding
# singular vectors of a sparse or dense matrix using a method of Baglama and
# Reichel.  It is a fast and memory-efficient way to compute a partial SVD.
#
# @param A numeric real- or complex-valued matrix or real-valued sparse matrix.
# @param nv number of right singular vectors to estimate.
# @param nu number of left singular vectors to estimate (defaults to \code{nv}).
# @param maxit maximum number of iterations.
# @param work working subspace dimension, larger values can speed convergence at the cost of more memory use.
# @param reorth if \code{TRUE}, apply full reorthogonalization to both SVD bases, otherwise
#   only apply reorthogonalization to the right SVD basis vectors; the latter case is cheaper per
#   iteration but, overall, may require more iterations for convergence. Automatically \code{TRUE}
#   when \code{fastpath=TRUE} (see below).
# @param tol convergence is determined when \eqn{\|A^TU - VS\| < tol\|A\|}{||A^T U - VS|| < tol*||A||},
#   and when the maximum relative change in estimated singular values from one iteration to the
#   next is less than \code{svtol = tol} (see \code{svtol} below),
#   where the spectral norm ||A|| is approximated by the
#   largest estimated singular value, and U, V, S are the matrices corresponding
#   to the estimated left and right singular vectors, and diagonal matrix of
#   estimated singular values, respectively.
# @param v optional starting vector or output from a previous run of \code{irlba} used
#   to restart the algorithm from where it left off (see the notes).
# @param right_only logical value indicating return only the right singular vectors
#  (\code{TRUE}) or both sets of vectors (\code{FALSE}). The right_only option can be
#  cheaper to compute and use much less memory when \code{nrow(A) >> ncol(A)} but note
#  that obtained solutions typically lose accuracy due to lack of re-orthogonalization in the
#  algorithm and that \code{right_only = TRUE} sets \code{fastpath = FALSE} (only use this option
#  for really large problems that run out of memory and when \code{nrow(A) >> ncol(A)}).
#  Consider increasing the \code{work} option to improve accuracy with \code{right_only=TRUE}.
# @param verbose logical value that when \code{TRUE} prints status messages during the computation.
# @param scale optional column scaling vector whose values divide each column of \code{A};
#   must be as long as the number of columns of \code{A} (see notes).
# @param center optional column centering vector whose values are subtracted from each
#   column of \code{A}; must be as long as the number of columns of \code{A} and may
#   not be used together with the deflation options below (see notes).
# @param shift optional shift value (square matrices only, see notes).
# @param mult DEPRECATED optional custom matrix multiplication function (default is \code{\%*\%}, see notes).
# @param fastpath try a fast C algorithm implementation if possible; set \code{fastpath=FALSE} to use the
#     reference R implementation. See the notes for more details.
# @param svtol additional stopping tolerance on maximum allowed absolute relative change across each
# estimated singular value between iterations.
# The default value of this parameter is to set it to \code{tol}. You can set \code{svtol=Inf} to
# effectively disable this stopping criterion. Setting \code{svtol=Inf} allows the method to
# terminate on the first Lanczos iteration if it finds an invariant subspace, but with less certainty
# that the converged subspace is the desired one. (It may, for instance, miss some of the largest
# singular values in difficult problems.)
# @param smallest set \code{smallest=TRUE} to estimate the smallest singular values and associated
# singular vectors. WARNING: this option is somewhat experimental, and may produce poor
# estimates for ill-conditioned matrices.
# @param ... optional additional arguments used to support experimental and deprecated features.
#
# @return
# Returns a list with entries:
# \describe{
#   \item{d:}{ max(nu, nv) approximate singular values}
#   \item{u:}{ nu approximate left singular vectors (only when right_only=FALSE)}
#   \item{v:}{ nv approximate right singular vectors}
#   \item{iter:}{ The number of Lanczos iterations carried out}
#   \item{mprod:}{ The total number of matrix vector products carried out}
# }
#
# @note
# The syntax of \code{irlba} partially follows \code{svd}, with an important
# exception. The usual R \code{svd} function always returns a complete set of
# singular values, even if the number of singular vectors \code{nu} or \code{nv}
# is set less than the maximum. The \code{irlba} function returns a number of
# estimated singular values equal to the maximum of the number of specified
# singular vectors \code{nu} and \code{nv}.
#
# Use the optional \code{scale} parameter to implicitly scale each column of
# the matrix \code{A} by the values in the \code{scale} vector, computing the
# truncated SVD of the column-scaled \code{sweep(A, 2, scale, FUN=`/`)}, or
# equivalently, \code{A \%*\% diag(1 / scale)}, without explicitly forming the
# scaled matrix. \code{scale} must be a non-zero vector of length equal
# to the number of columns of \code{A}.
#
# Use the optional \code{center} parameter to implicitly subtract the values
# in the \code{center} vector from each column of \code{A}, computing the
# truncated SVD of \code{sweep(A, 2, center, FUN=`-`)},
# without explicitly forming the centered matrix. \code{center}
# must be a vector of length equal to the number of columns of \code{A}.
# This option may be used to efficiently compute principal components without
# explicitly forming the centered matrix (which can, importantly, preserve
# sparsity in the matrix). See the examples.
#
# The optional \code{shift} scalar valued argument applies only to square matrices; use it
# to estimate the partial svd of \code{A + diag(shift, nrow(A), nrow(A))}
# (without explicitly forming the shifted matrix).
#
# (Deprecated) Specify an optional alternative matrix multiplication operator in the
# \code{mult} parameter. \code{mult} must be a function of two arguments,
# and must handle both cases where one argument is a vector and the other
# a matrix. This option is deprecated and will be removed in a future version.
# The new preferred method simply uses R itself to define a custom matrix class
# with your user-defined matrix multiplication operator. See the examples.
#
# Use the \code{v} option to supply a starting vector for the iterative
# method. A random vector is used by default (precede with \code{set.seed()}
# for reproducibility). Optionally set \code{v} to
# the output of a previous run of \code{irlba} to restart the method, adding
# additional singular values/vectors without recomputing the solution
# subspace. See the examples.
#
# The function may generate the following warnings:
# \itemize{
#   \item{"did not converge--results might be invalid!; try increasing work or maxit"
#   means that the algorithm didn't
#   converge -- this is potentially a serious problem and the returned results may not be valid. \code{irlba}
#   reports a warning here instead of an error so that you can inspect whatever is returned. If this
#   happens, carefully heed the warning and inspect the result. You may also try setting \code{fastpath=FALSE}.}
#   \item{"You're computing a large percentage of total singular values, standard svd might work better!"
#     \code{irlba} is designed to efficiently compute a few of the largest singular values and associated
#      singular vectors of a matrix. The standard \code{svd} function will be more efficient for computing
#      large numbers of singular values than \code{irlba}.}
#    \item{"convergence criterion below machine epsilon" means that the product of \code{tol} and the
#      largest estimated singular value is really small and the normal convergence criterion is only
#      met up to round off error.}
# }
# The function might return an error for several reasons including a situation when the starting
# vector \code{v} is near the null space of the matrix. In that case, try a different \code{v}.
#
# The \code{fastpath=TRUE} option only supports real-valued matrices and sparse matrices
# of type \code{dgCMatrix} (for now). Other problems fall back to the reference
# R implementation.
#
# @references
# Baglama, James, and Lothar Reichel. "Augmented implicitly restarted Lanczos bidiagonalization methods." SIAM Journal on Scientific Computing 27.1 (2005): 19-42.
#
# @examples
# set.seed(1)
#
# A <- matrix(runif(400), nrow=20)
# S <- irlba(A, 3)
# S$d
#
# # Compare with svd
# svd(A)$d[1:3]
#
# # Restart the algorithm to compute more singular values
# # (starting with an existing solution S)
# S1 <- irlba(A, 5, v=S)
#
# # Estimate smallest singular values
# irlba(A, 3, smallest=TRUE)$d
#
# #Compare with
# tail(svd(A)$d, 3)
#
# # Principal components (see also prcomp_irlba)
# P <- irlba(A, nv=1, center=colMeans(A))
#
# # Compare with prcomp and prcomp_irlba (might vary up to sign)
# cbind(P$v,
#       prcomp(A)$rotation[, 1],
#       prcomp_irlba(A)$rotation[, 1])
#
# # A custom matrix multiplication function that scales the columns of A
# # (cf the scale option). This function scales the columns of A to unit norm.
# col_scale <- sqrt(apply(A, 2, crossprod))
# setClass("scaled_matrix", contains="matrix", slots=c(scale="numeric"))
# setMethod("%*%", signature(x="scaled_matrix", y="numeric"),
#    function(x ,y) x@.Data %*% (y / x@scale))
# setMethod("%*%", signature(x="numeric", y="scaled_matrix"),
#    function(x ,y) (x %*% y@.Data) / y@scale)
# a <- new("scaled_matrix", A, scale=col_scale)
# irlba(a, 3)$d
#
# # Compare with:
# svd(sweep(A, 2, col_scale, FUN=`/`))$d[1:3]
#
#
# @seealso \code{\link{svd}}, \code{\link{prcomp}}, \code{\link{partial_eigen}}, \code{\link{svdr}}
# @import Matrix
# @importFrom stats rnorm
# @importFrom methods slotNames
irlba__irlba <-
function(A,                     # data matrix
         nv=5, nu=nv,           # number of singular vectors to estimate
         maxit=1000,            # maximum number of iterations
         work=nv + 7,           # working subspace size
         reorth=TRUE,           # TRUE=full reorthogonalization
         tol=1e-5,              # stopping tolerance
         v=NULL,                # optional starting vector or restart
         right_only=FALSE,      # TRUE=only return V
         verbose=FALSE,         # display status messages
         scale=NULL,            # optional column scaling
         center=NULL,           # optional column centering
         shift=NULL,            # optional shift for square matrices
         mult=NULL,             # optional custom matrix multiplication func.
         fastpath=TRUE,         # use the faster C implementation if possible
         svtol=tol,             # stopping tolerance percent change in estimated svs
         smallest=FALSE,        # set to TRUE to estimate subspaces associated w/smallest singular values
         ...)                   # optional experimental or deprecated arguments
{
# ---------------------------------------------------------------------
# Check input parameters
# ---------------------------------------------------------------------
  ropts <- options(warn=1) # immediately show warnings
  mflag <- new.env()
  mflag$flag <- FALSE
  on.exit(options(ropts))
  interchange <- FALSE
  eps <- .Machine$double.eps
  # hidden support for old, removed (previously deprecated) parameters
  # this is here as a convenience to keep old code working without change
  # also supports experimental features not yet promoted to the api
  mcall <- as.list(match.call())
  random <- eval(mcall[["rng"]])
  if (is.null(random)) random <- stats::rnorm # default RNG
  # Maximum number of Ritz vectors to use in augmentation, may be less
  # depending on workspace size.
  maxritz <- eval(mcall[["maxritz"]]) # experimental
  if (is.null(maxritz)) maxritz <- 3
  eps2 <- eval(mcall[["invariant_subspace_tolerance"]])
  if (is.null(eps2)) eps2 <- eps ^ (4 / 5)
  du <- eval(mcall[["du"]]) # deprecated
  dv <- eval(mcall[["dv"]]) # deprecated
  ds <- eval(mcall[["ds"]]) # deprecated
  deflate <- is.null(du) + is.null(ds) + is.null(dv)
  if (is.logical(scale) && ! scale) scale <- NULL
  if (is.logical(shift) && ! shift) shift <- NULL
  if (is.logical(center) && ! center) center <- NULL
  if (smallest) fastpath <- FALSE  # for now anyway
  if (any(dim(A) > 2 ^ 32 - 1)) fastpath <- FALSE # for now
  if (deflate == 3)
  {
    deflate <- FALSE
  } else if (deflate == 0)
  {
    deflate <- TRUE
    warning("The deflation options have been deprecated. Please modify your code to not use them.")
    if (length(ds) > 1) stop("deflation limited to one dimension")
    if (!is.null(dim(du))) du <- du[, 1]
    if (!is.null(dim(dv))) dv <- dv[, 1]
  } else stop("all three du ds dv parameters must be specified for deflation")
  if (!is.null(center))
  {
    if (is.logical(center) && center) center <- colMeans(A)
    if (deflate) stop("the center parameter can't be specified together with deflation parameters")
    if (length(center) != ncol(A)) stop("center must be a vector of length ncol(A)")
    if (fastpath && ! right_only) du <- NULL
    else du <- 1
    ds <- 1
    dv <- center
    deflate <- TRUE
  }
  if ("integer" == typeof(A)) A <- A + 0.0
  iscomplex <- is.complex(A)
  m <- nrow(A)
  n <- ncol(A)
  if (is.null(nu)) nu <- nv
  if (!is.null(mult) && deflate) stop("the mult parameter can't be specified together with deflation parameters")
  missingmult <- FALSE
  if (is.null(mult))
  {
    missingmult <- TRUE
    mult <- `%*%`
  }
  k <- max(nu, nv)
  if (k <= 0)  stop("max(nu, nv) must be positive")
  if (k > min(m - 1, n - 1)) stop("max(nu, nv) must be strictly less than min(nrow(A), ncol(A))")
  if (k >= 0.5 * min(m, n))
  {
    warning("You're computing too large a percentage of total singular values, use a standard svd instead.")
  }
  if (work <= 1) stop("work must be greater than 1")
  if (tol < 0) stop("tol must be non-negative")
  if (maxit <= 0) stop("maxit must be positive")
  # work must be strictly larger than requested subspace dimension, except see right_only below
  if (work <= k && ! right_only) work <- k + 1
  if (work >= min(n, m))
  {
    work <- min(n, m)
    if (work <= k)
    {
      k <- work - 1  # the best we can do! Need to reduce output subspace dimension
      warning("Requested subspace dimension too large! Reduced to ", k)
    }
  }
  k_org <- k
  w_dim <- work
  if (right_only)
  {
    w_dim <- 1
    fastpath <- FALSE
  }
  if (n > m && smallest)
  {
    # Interchange dimensions m,n so that dim(A'A) = min(m,n) when seeking the
    # smallest singular values; avoids finding zero-valued smallest singular values.
    interchange <- TRUE
    temp <- m
    m <- n
    n <- temp
  }

  if (verbose)
  {
    message("Working dimension size ", work)
  }
# Check for tiny problem, use standard SVD in that case. Make definition of 'tiny' larger?
  if (min(m, n) < 6)
  {
    A <- as.matrix(A) # avoid need to define "+" and "/" for arbitrary matrix types.
    if (verbose) message("Tiny problem detected, using standard `svd` function.")
    if (!is.null(scale)) {
      A <- sweep(A, 2, scale, "/")
      dv <- dv / scale # scale the centering vector.
    }
    if (!is.null(shift)) A <- A + diag(shift, nrow(A), ncol(A))
    if (deflate)
    {
      if (is.null(du)) du <- rep(1, nrow(A))
      A <- A - (ds * du) %*% t(dv)
    }
    s <- svd(A)
    if (smallest)
    {
      return(list(d=tail(s$d, k), u=s$u[, tail(seq(ncol(s$u)), k), drop=FALSE],
              v=s$v[, tail(seq(ncol(s$v), k)), drop=FALSE], iter=0, mprod=0))
    }
    return(list(d=s$d[1:k], u=s$u[, 1:nu, drop=FALSE],
              v=s$v[, 1:nv, drop=FALSE], iter=0, mprod=0))
  }

# Try to use the fast C-language code path
  if (deflate) fastpath <- fastpath && is.null(du)
# Only matrix, dgCMatrix supported by fastpath
  fastpath <- fastpath && (("Matrix" %in% attributes(class(A)) && ("dgCMatrix" %in% class(A))) || "matrix" %in% class(A))

# Allocate memory for W and F:
  W <- matrix(0.0, m, w_dim)
  F <- matrix(0.0, n, 1)
  restart <- FALSE
  if (is.list(v))
  {
    if (is.null(v$v) || is.null(v$d) || is.null(v$u)) stop("restart requires left and right singular vectors")
    if (max(nu, nv) <= min(ncol(v$u), ncol(v$v))) return(v) # Nothing to do!
    right_only <- FALSE
    W[, 1:ncol(v$u)] <- v$u
    d <- v$d
    V <- matrix(0.0, n, work)
    V[, 1:ncol(v$v)] <- v$v
    restart <- TRUE
  } else if (is.null(v))
  {
# If starting matrix v is not given then set V to be an (n x 1) matrix of
# normally distributed random numbers.  In any case, allocate V appropriate to
# problem size:
    V <- matrix(0.0, n, work)
    V[, 1] <- random(n)
  } else
  {
# user-supplied starting subspace
    V <- matrix(0.0, n, work)
    V[1:length(v)] <- v
  }

# ---------------------------------------------------------------------
# Initialize local variables
# ---------------------------------------------------------------------
  B <- NULL                  # Bidiagonal matrix
  Bsz <- NULL                # Size of B
  eps23 <- eps ^ (2 / 3)     # Used for Smax/avoids using zero
  iter <- 1                  # Man loop iteration count
  mprod <- 0                 # Number of matrix-vector products
  R_F <- NULL                # 2-norm of residual vector F
  sqrteps <- sqrt(eps)       #
  Smax <- 1                  # Max value of all computed singular values of
                             # B est. ||A||_2
  Smin <- NULL               # Min value of all computed singular values of
                             # B est. cond(A)
  lastsv <- c()              # estimated sv in last iteration

# Check for user-supplied restart condition
  if (restart)
  {
    B <- cbind(diag(d), 0)
    k <- length(d)

    F <- random(n)
    F <- orthog(F, V[, 1:k])
    V[, k + 1] <- F / norm2(F)
  }

# Change du to be non-NULL, for non-fastpath'able matrices with non-NULL scale.
  if (deflate && is.null(du)) du <- 1

# ---------------------------------------------------------------------
# Main iteration
# ---------------------------------------------------------------------
  while (iter <= maxit)
  {
# ---------------------------------------------------------------------
# Compute the Lanczos bidiagonal decomposition:
# such that AV  = WB
# and       t(A)W = VB + Ft(E)
# This routine updates W, V, F, B, mprod
#
# Note on scale and center: These options are applied implicitly below
# for maximum computational efficiency. This complicates their application
# somewhat, but saves a few flops.
# ---------------------------------------------------------------------
    j <- 1
#   Normalize starting vector:
    if (iter == 1 && !restart)
    {
      V[, 1] <- V[, 1] / norm2(V[, 1])
    }
    else j <- k + 1
#   j_w is used here to support the right_only=TRUE case.
    j_w <- ifelse(w_dim > 1, j, 1)

#   Compute W=AV
#   Optionally apply scale
    VJ <- V[, j]
    if (!is.null(scale))
    {
      VJ <- VJ / scale
    }
    if (interchange) avj <- mult(VJ, A)
    else avj <- mult(A, VJ)

#   Handle non-ordinary arrays as products.
    W[, j_w] <- as.vector(avj)
    mprod <- mprod + 1

#   Optionally apply shift
    if (!is.null(shift))
    {
      W[, j_w] <- W[, j_w] + shift * VJ
    }

#   Optionally apply deflation
    if (deflate)
    {
      W[, j_w] <- W[, j_w] - ds * drop(cross(dv, VJ)) * du
    }

#   Orthogonalize W
    if (iter != 1 && w_dim > 1 && reorth)
    {
      W[, j] <- orthog(W[, j, drop=FALSE], W[, 1:(j - 1), drop=FALSE])
    }

    S <- norm2(W[, j_w, drop=FALSE])
#   Check for linearly dependent vectors
    if (is.na(S) || S < eps2 && j == 1) stop("starting vector near the null space")
    if (is.na(S) || S < eps2)
    {
      if (verbose) message_once("invariant subspace found", flag=mflag)
      W[, j_w] <- random(nrow(W))
      if (w_dim > 1) W[, j] <- orthog(W[, j], W[, 1:(j - 1)])
      W[, j_w] <- W[, j_w] / norm2(W[, j_w])
      S <- 0
    }
    else W[, j_w] <- W[, j_w] / S

#   Lanczos process
    while (j <= work)
    {
      j_w <- ifelse(w_dim > 1, j, 1)
      if (iscomplex)
      {
        if (interchange) F <- Conj(t(drop(mult(A, Conj(drop(W[, j_w]))))))
        else F <- Conj(t(drop(mult(Conj(drop(W[, j_w])), A))))
      }
      else
      {
        if (interchange) F <- t(drop(mult(A, drop(W[, j_w]))))
        else F <- t(drop(mult(drop(W[, j_w]), A)))
      }
#     Optionally apply shift, scale, deflate
      if (!is.null(shift)) F <- F + shift * W[, j_w]
      if (!is.null(scale)) F <- F / scale
      if (deflate) {
        sub <- sum(W[, j_w]) * dv
        if (!is.null(scale)) sub <- sub / scale
        F <- F - sub
      }
      mprod <- mprod + 1
      F <- drop(F - S * V[, j])
#     Orthogonalize
      F <- orthog(F, V[, 1:j, drop=FALSE])
      if (j + 1 <= work)
      {
        R <- norm2(F)
#       Check for linear dependence
        if (R < eps2)
        {
          if (verbose) message_once("invariant subspace found", flag=mflag)
          F <- matrix(random(dim(V)[1]), dim(V)[1], 1)
          F <- orthog(F, V[, 1:j, drop=FALSE])
          V[, j + 1] <- F / norm2(F)
          R <- 0
        }
        else V[, j + 1] <- F / R

#       Compute block diagonal matrix
        if (is.null(B)) B <- cbind(S, R)
        else            B <- rbind(cbind(B, 0), c(rep(0, ncol(B) - 1), S, R))

        jp1_w <- ifelse(w_dim > 1, j + 1, 1)
        w_old <- W[, j_w]

#       Optionally apply scale
        VJP1 <- V[, j + 1]
        if (!is.null(scale))
        {
          VJP1 <- VJP1 / scale
        }
        if (interchange) W[, jp1_w] <- drop(mult(drop(VJP1), A))
        else W[, jp1_w] <- drop(mult(A, drop(VJP1)))
        mprod <- mprod + 1

#       Optionally apply shift
        if (!is.null(shift))
        {
          W[, jp1_w] <- W[, jp1_w] + shift * VJP1
        }

#       Optionally apply deflation
        if (deflate)
        {
          W[, jp1_w] <- W[, jp1_w] - ds * drop(cross(dv, VJP1)) * du
        }

#       One step of the classical Gram-Schmidt process
        W[, jp1_w] <- W[, jp1_w] - R * w_old

#       Full reorthogonalization of W
        if (reorth && w_dim > 1) W[, j + 1] <- orthog(W[, j + 1], W[, 1:j])
        S <- norm2(W[, jp1_w])
#       Check for linear dependence
        if (S < eps2)
        {
          if (verbose) message_once("invariant subspace found", flag=mflag)
          W[, jp1_w] <- random(nrow(W))
          if (w_dim > 1) W[, j + 1] <- orthog(W[, j + 1], W[, 1:j])
          W[, jp1_w] <- W[, jp1_w] / norm2(W[, jp1_w])
          S <- 0
        }
        else W[, jp1_w] <- W[, jp1_w] / S
      }
      else
      {
#       Add a last block to matrix B
        B <- rbind(B, c(rep(0, j - 1), S))
      }
      j <- j + 1
    }
# ---------------------------------------------------------------------
# (End of the Lanczos bidiagonalization part)
# ---------------------------------------------------------------------
    Bsz <- nrow(B)
    R_F <- norm2(F)
    F <- F / R_F
#   Compute singular triplets of B, svd must return ordered singular
#   values from largest to smallest.
    Bsvd <- svd(B)

#   Estimate ||A|| using the largest singular value over all iterations
#   and estimate the cond(A) using approximations to the largest and
#   smallest singular values. If a small singular value is less than sqrteps
#   require two-sided reorthogonalization.
    if (iter == 1)
    {
      Smax <- Bsvd$d[1]
      Smin <- Bsvd$d[Bsz]
    }
    else
    {
      Smax <- max(Smax, Bsvd$d[1])
      Smin <- min(Smin, Bsvd$d[Bsz])
    }
    Smax <- max(eps23, Smax)
    if (! reorth && Smin / Smax < sqrteps)
    {
      warning("The matrix is ill-conditioned. Basis will be reorthogonalized.")
      reorth <- TRUE
    }
    if (smallest)
    {
      jj <- seq(ncol(Bsvd$u), 1, by = -1)
      Bsvd$u <- Bsvd$u[, jj]
      Bsvd$d <- Bsvd$d[jj]
      Bsvd$v <- Bsvd$v[, jj]
    }

#   Compute the residuals
    R <- R_F * Bsvd$u[Bsz, , drop=FALSE]
#   Check for convergence
    ct <- convtests(Bsz, tol, k_org, Bsvd, abs(R), k, Smax, lastsv, svtol, maxritz, work, S)
    if (verbose)
    {
      message("iter= ", iter,
              ", mprod= ", mprod,
              ", sv[", k_org, "]=", sprintf("%.2e", Bsvd$d[k_org]),
              ", %change=", sprintf("%.2e", (Bsvd$d[k_org] - lastsv[k_org])/Bsvd$d[k_org]),
              ", k=", ct$k)
    }
    lastsv <- Bsvd$d
    k <- ct$k

#   If all desired singular values converged, then exit main loop
    if (ct$converged) break
    if (iter >= maxit) break

#   Compute the starting vectors and first block of B
    if (smallest && (Smin / Smax > sqrteps))
    {
#     Update the SVD of B to be the svd of [B ||F||E_m]
      Bsvd2.d <- Bsvd$d
      Bsvd2.d <- diag(Bsvd2.d, nrow=length(Bsvd2.d))
      Bsvd2 <- svd(cbind(Bsvd2.d, t(R)))
      jj <- seq(ncol(Bsvd2$u), 1, by=-1)
      Bsvd2$u <- Bsvd2$u[, jj]
      Bsvd2$d <- Bsvd2$d[jj]
      Bsvd2$v <- Bsvd2$v[, jj]

      Bsvd$d <- Bsvd2$d
      Bsvd$u <- Bsvd$u %*% Bsvd2$u
      Bsvd$v <- cbind(rbind(Bsvd$v, rep(0, Bsz)), c(rep(0, Bsz), 1)) %*% Bsvd2$v
      V_B_last <- Bsvd$v[Bsz + 1, 1:k]
      s <- R_F * solve(B, cbind(c(rep(0, Bsz - 1), 1)))
      Bsvd$v <- Bsvd$v[1:Bsz, , drop=FALSE] + s %*% Bsvd$v[Bsz + 1, ]

      qrv <- qr(cbind(rbind(Bsvd$v[, 1:k], 0), rbind(-s, 1)))
      Bsvd$v <- qr.Q(qrv)
      R <- qr.R(qrv)
      V[, 1:(k + 1)] <- cbind(V, F) %*% Bsvd$v

#  Update and compute the k by k+1 part of B
      UT <- t(R[1:(k + 1), 1:k] + R[, k + 1] %*% rbind(V_B_last))
      B <- diag(Bsvd$d[1:k], nrow=k) %*% (UT * upper.tri(UT, diag=TRUE))[1:k, 1:(k+1)]
    } else
    {
#   using the Ritz vectors
      V[, 1:(k + dim(F)[2])] <- cbind(V[, 1:(dim(Bsvd$v)[1]), drop=FALSE] %*% Bsvd$v[, 1:k], F)
      B <- cbind(diag(Bsvd$d[1:k], nrow=k), R[1:k])
    }

#   Update the left approximate singular vectors
    if (w_dim > 1)
    {
      W[, 1:k] <- W[, 1:(dim(Bsvd$u)[1]), drop=FALSE] %*% Bsvd$u[, 1:k]
    }

    iter <- iter + 1
  }
# ---------------------------------------------------------------------
# End of the main iteration loop
# Output results
# ---------------------------------------------------------------------
  if (!ct$converged) warning("did not converge--results might be invalid!; try increasing maxit or work")
  d <- Bsvd$d[1:k_org]
  if (!right_only)
  {
    u <- W[, 1:(dim(Bsvd$u)[1]), drop=FALSE] %*% Bsvd$u[, 1:k_org, drop=FALSE]
  }
  v <- V[, 1:(dim(Bsvd$v)[1]), drop=FALSE] %*% Bsvd$v[, 1:k_org, drop=FALSE]
  if (smallest)
  {
    reverse <- seq(length(d), 1)
    d <- d[reverse]
    if (!right_only) u <- u[, reverse]
    v <- v[, reverse]
  }
  if (tol * d[1] < eps) warning("convergence criterion below machine epsilon")
  if (right_only)
    return(list(d=d, v=v[, 1:nv, drop=FALSE], iter=iter, mprod=mprod))
  return(list(d=d, u=u[, 1:nu, drop=FALSE],
              v=v[, 1:nv, drop=FALSE], iter=iter, mprod=mprod))
}

# Find a few approximate largest eigenvalues and corresponding eigenvectors of a symmetric matrix.
#
# Use \code{partial_eigen} to estimate a subset of the largest (most positive)
# eigenvalues and corresponding eigenvectors of a symmetric dense or sparse
# real-valued matrix.
#
# @param x numeric real-valued dense or sparse matrix.
# @param n number of largest eigenvalues and corresponding eigenvectors to compute.
# @param symmetric \code{TRUE} indicates \code{x} is a symmetric matrix (the default);
#   specify \code{symmetric=FALSE} to compute the largest eigenvalues and corresponding
#   eigenvectors of \code{t(x) \%*\% x} instead.
# @param ... optional additional parameters passed to the \code{irlba} function.
#
# @return
# Returns a list with entries:
# \itemize{
#   \item{values}{ n approximate largest eigenvalues}
#   \item{vectors}{ n approximate corresponding eigenvectors}
# }
#
# @note
# Specify \code{symmetric=FALSE} to compute the largest \code{n} eigenvalues
# and corresponding eigenvectors of the symmetric matrix cross-product
# \code{t(x) \%*\% x}.
#
# This function uses the \code{irlba} function under the hood. See \code{?irlba}
# for description of additional options, especially the \code{tol} parameter.
#
# See the RSpectra package https://cran.r-project.org/package=RSpectra for more comprehensive
# partial eigenvalue decomposition.
#
# @references
# Augmented Implicitly Restarted Lanczos Bidiagonalization Methods, J. Baglama and L. Reichel, SIAM J. Sci. Comput. 2005.
#
# @examples
# set.seed(1)
# # Construct a symmetric matrix with some positive and negative eigenvalues:
# V <- qr.Q(qr(matrix(runif(100), nrow=10)))
# x <- V %*% diag(c(10, -9, 8, -7, 6, -5, 4, -3, 2, -1)) %*% t(V)
# partial_eigen(x, 3)$values
#
# # Compare with eigen
# eigen(x)$values[1:3]
#
# # Use symmetric=FALSE to compute the eigenvalues of t(x) %*% x for general
# # matrices x:
# x <- matrix(rnorm(100), 10)
# partial_eigen(x, 3, symmetric=FALSE)$values
# eigen(crossprod(x))$values
#
# @seealso \code{\link{eigen}}, \code{\link{irlba}}
irlba__partial_eigen <- function(x, n=5, symmetric=TRUE, ...)
{
  if (n > 0.5 * min(nrow(x), ncol(x)))
  {
    warning("You're computing a large percentage of total eigenvalues, the standard eigen function will likely work better!")
  }
  if (!symmetric)
  {
    L <- irlba__irlba(x, n, ...)
    return(list(vectors=L$v, values=L$d ^ 2))
  }
  L <- irlba__irlba(x, n, ...)
  s <- sign(L$u[1, ] * L$v[1, ])
  if (all(s > 0))
  {
    return(list(vectors=L$u, values=L$d))
  }
  i <- min(which(s < 0))
  shift <- L$d[i]
  L <- irlba__irlba(x, n, shift=shift, ...)
  return(list(vectors=L$u, values=L$d - shift))
}

# Principal Components Analysis
#
# Efficient computation of a truncated principal components analysis of a given data matrix
# using an implicitly restarted Lanczos method from the \code{\link{irlba}} package.
#
# @param x a numeric or complex matrix (or data frame) which provides
#          the data for the principal components analysis.
# @param retx a logical value indicating whether the rotated variables should be returned.
# @param center a logical value indicating whether the variables should be
#          shifted to be zero centered. Alternately, a centering vector of length
#          equal the number of columns of \code{x} can be supplied.
# @param scale. a logical value indicating whether the variables should be
#          scaled to have unit variance before the analysis takes place.
#          The default is \code{FALSE} for consistency with S, but scaling is often advisable.
#          Alternatively, a vector of length equal the number of columns of \code{x} can be supplied.
#
#          The value of \code{scale} determines how column scaling is performed
#          (after centering).  If \code{scale} is a numeric vector with length
#          equal to the number of columns of \code{x}, then each column of \code{x} is
#          divided by the corresponding value from \code{scale}.  If \code{scale} is
#          \code{TRUE} then scaling is done by dividing the (centered) columns of
#          \code{x} by their standard deviations if \code{center=TRUE}, and the
#          root mean square otherwise.  If \code{scale} is \code{FALSE}, no scaling is done.
#          See \code{\link{scale}} for more details.
# @param n integer number of principal component vectors to return, must be less than
# \code{min(dim(x))}.
# @param ... additional arguments passed to \code{\link{irlba}}.
#
# @return
# A list with class "prcomp" containing the following components:
# \itemize{
#    \item{sdev} {the standard deviations of the principal components (i.e.,
#          the square roots of the eigenvalues of the
#          covariance/correlation matrix, though the calculation is
#          actually done with the singular values of the data matrix).}
#   \item{rotation} {the matrix of variable loadings (i.e., a matrix whose columns
#          contain the eigenvectors).}
#   \item {x} {if \code{retx} is \code{TRUE} the value of the rotated data (the centred
#          (and scaled if requested) data multiplied by the \code{rotation}
#         matrix) is returned.  Hence, \code{cov(x)} is the diagonal matrix
#          \code{diag(sdev^2)}.}
#   \item{center, scale} {the centering and scaling used, or \code{FALSE}.}
# }
#
# @note
# The signs of the columns of the rotation matrix are arbitrary, and
# so may differ between different programs for PCA, and even between
# different builds of R.
#
# NOTE DIFFERENCES WITH THE DEFAULT \code{\link{prcomp}} FUNCTION!
# The \code{tol} truncation argument found in \code{prcomp} is not supported.
# In place of the truncation tolerance in the original function, the
# \code{prcomp_irlba}  function has the argument \code{n} explicitly giving the
# number of principal components to return. A warning is generated if the
# argument \code{tol} is used, which is interpreted differently between
# the two functions.
#
# @examples
# set.seed(1)
# x  <- matrix(rnorm(200), nrow=20)
# p1 <- prcomp_irlba(x, n=3)
# summary(p1)
#
# # Compare with
# p2 <- prcomp(x, tol=0.7)
# summary(p2)
#
#
# @seealso \code{\link{prcomp}}
# @import Matrix
# @importFrom stats rnorm prcomp sd var
# @importFrom methods slotNames slot
irlba__prcomp_irlba <- function(x, n = 3, retx = TRUE, center = TRUE, scale. = FALSE, ...)
{
  a <- names(as.list(match.call()))
  ans <- list(scale=scale.)
  if ("tol" %in% a)
    warning("The `tol` truncation argument from `prcomp` is not supported by
`prcomp_irlba`. If specified, `tol` is passed to the `irlba` function to
control that algorithm's convergence tolerance. See `?prcomp_irlba` for help.")
  # Try to convert data frame to matrix...
  if (is.data.frame(x)) x <- as.matrix(x)
  args <- list(A=x, nv=n)
  if (is.logical(center))
  {
    if (center) args$center <- colMeans(x)
  } else args$center <- center
  if (is.logical(scale.))
  {
    if (is.numeric(args$center))
    {
      f <- function(i) sqrt(sum((x[, i] - args$center[i]) ^ 2) / (nrow(x) - 1L))
      scale. <- vapply(seq(ncol(x)), f, pi, USE.NAMES=FALSE)
      if (ans$scale) ans$totalvar <- ncol(x)
      else ans$totalvar <- sum(scale. ^ 2)
    } else
    {
      if (ans$scale)
      {
        scale. <- apply(x, 2L, function(v) sqrt(sum(v ^ 2) / max(1, length(v) - 1L)))
        f <- function(i) sqrt(sum((x[, i] / scale.[i]) ^ 2) / (nrow(x) - 1L))
        ans$totalvar <- sum(vapply(seq(ncol(x)), f, pi, USE.NAMES=FALSE) ^ 2)
      } else
      {
        f <- function(i) sum(x[, i] ^ 2) / (nrow(x) - 1L)
        ans$totalvar <- sum(vapply(seq(ncol(x)), f, pi, USE.NAMES=FALSE))
      }
    }
    if (ans$scale) args$scale <- scale.
  } else
  {
    args$scale <- scale.
    f <- function(i) sqrt(sum((x[, i] / scale.[i]) ^ 2) / (nrow(x) - 1L))
    ans$totalvar <- sum(vapply(seq(ncol(x)), f, pi, USE.NAMES=FALSE))
  }
  if (!missing(...)) args <- c(args, list(...))

  s <- do.call(irlba__irlba, args=args)
  ans$sdev <- s$d / sqrt(max(1, nrow(x) - 1))
  ans$rotation <- s$v
  colnames(ans$rotation) <- paste("PC", seq(1, ncol(ans$rotation)), sep="")
  ans$center <- args$center
  if (retx)
  {
    ans <- c(ans, list(x = sweep(s$u, 2, s$d, FUN=`*`)))
    colnames(ans$x) <- paste("PC", seq(1, ncol(ans$rotation)), sep="")
  }
  class(ans) <- c("irlba_prcomp", "prcomp")
  ans
}

# Summary method for truncated pca objects computed by \code{prcomp_irlba}.
# @param object An object returned by \code{prcomp_irlba}.
# @param ... Optional arguments passed to \code{summary}.
# @method summary irlba_prcomp
summary.irlba_prcomp <- function(object, ...)
{
  chkDots(...)
  vars <- object$sdev ^ 2
  vars <- vars / object$totalvar
  importance <- rbind("Standard deviation" = object$sdev,
                      "Proportion of Variance" = round(vars, 5),
                      "Cumulative Proportion" = round(cumsum(vars), 5))
  k <- ncol(object$rotation)
  colnames(importance) <- c(colnames(object$rotation), rep("", length(vars) - k))
  object$importance <- importance
  class(object) <- "summary.prcomp"
  object
}

# Sparse regularized low-rank matrix approximation.
#
# Estimate an \eqn{{\ell}1}{l1}-penalized
# singular value or principal components decomposition (SVD or PCA) that introduces sparsity in the
# right singular vectors based on the fast and memory-efficient
# sPCA-rSVD algorithm of Haipeng Shen and Jianhua Huang.
# @param x A numeric real- or complex-valued matrix or real-valued sparse matrix.
# @param k Matrix rank of the computed decomposition (see the Details section below).
# @param n Number of nonzero components in the right singular vectors. If \code{k > 1},
#        then a single value of \code{n} specifies the number of nonzero components
#        in each regularized right singular vector. Or, specify a vector of length
#        \code{k} indicating the number of desired nonzero components in each
#        returned vector. See the examples.
# @param maxit Maximum number of soft-thresholding iterations.
# @param tol Convergence is determined when \eqn{\|U_j - U_{j-1}\|_F < tol}{||U_j - U_{j-1}||_F < tol}, where \eqn{U_j} is the matrix of estimated left regularized singular vectors at iteration \eqn{j}.
# @param center a logical value indicating whether the variables should be
#        shifted to be zero centered. Alternately, a centering vector of length
#        equal the number of columns of \code{x} can be supplied. Use \code{center=TRUE}
#        to perform a regularized sparse PCA.
# @param scale. a logical value indicating whether the variables should be
#        scaled to have unit variance before the analysis takes place.
#        Alternatively, a vector of length equal the number of columns of \code{x} can be supplied.
#
#        The value of \code{scale} determines how column scaling is performed
#        (after centering).  If \code{scale} is a numeric vector with length
#        equal to the number of columns of \code{x}, then each column of \code{x} is
#        divided by the corresponding value from \code{scale}.  If \code{scale} is
#        \code{TRUE} then scaling is done by dividing the (centered) columns of
#        \code{x} by their standard deviations if \code{center=TRUE}, and the
#        root mean square otherwise.  If \code{scale} is \code{FALSE}, no scaling is done.
#        See \code{\link{scale}} for more details.
# @param alpha Optional  scalar regularization parameter between zero and one (see Details below).
# @param  tsvd Optional initial rank-k truncated SVD or PCA (skips computation if supplied).
# @param ... Additional arguments passed to \code{\link{irlba}}.
# @details
# The \code{ssvd} function implements a version of an algorithm by
# Shen and Huang that computes a penalized SVD or PCA that introduces
# sparsity in the right singular vectors by solving a penalized least squares problem.
# The algorithm in the rank 1 case finds vectors \eqn{u, w}{u, w} that minimize
# \deqn{\|x - u w^T\|_F^2 + \lambda \|w\|_1}{||x - u w^T||_F^2 + lambda||w||_1}
# such that \eqn{\|u\| = 1}{||u|| = 1},
# and then sets \eqn{v = w / \|w\|}{v = w / ||w||} and
# \eqn{d = u^T x v}{d = u^T x v};
# see the referenced paper for details. The penalty \eqn{\lambda}{lambda} is
# implicitly determined from the specified desired number of nonzero values \code{n}.
# Higher rank output is determined similarly
# but using a sequence of \eqn{\lambda}{lambda} values determined to maintain the desired number
# of nonzero elements in each column of \code{v} specified by \code{n}.
# Unlike standard SVD or PCA, the columns of the returned \code{v} when \code{k > 1} may not be orthogonal.
#
# @return
# A list containing the following components:
# \itemize{
#    \item{u} {regularized left singular vectors with orthonormal columns}
#    \item{d} {regularized upper-triangluar projection matrix so that \code{x \%*\% v == u \%*\% d}}
#    \item{v} {regularized, sparse right singular vectors with columns of unit norm}
#    \item{center, scale} {the centering and scaling used, if any}
#    \item{lambda} {the per-column regularization parameter found to obtain the desired sparsity}
#    \item{iter} {number of soft thresholding iterations}
#    \item{n} {value of input parameter \code{n}}
#    \item{alpha} {value of input parameter \code{alpha}}
# }
# @note
# Our \code{ssvd} implementation of the Shen-Huang method makes the following choices:
# \enumerate{
# \item{The l1 penalty is the only available penalty function. Other penalties may appear in the future.}
# \item{Given a desired number of nonzero elements in \code{v}, value(s) for the \eqn{\lambda}{lambda}
#       penalty are determined to achieve the sparsity goal subject to the parameter \code{alpha}.}
# \item{An experimental block implementation is used for results with rank greater than 1 (when \code{k > 1})
#       instead of the deflation method described in the reference.}
# \item{The choice of a penalty lambda associated with a given number of desired nonzero
#       components is not unique. The \code{alpha} parameter, a scalar between zero and one,
#       selects any possible value of lambda that produces the desired number of
#       nonzero entries. The default \code{alpha = 0} selects a penalized solution with
#       largest corresponding value of \code{d} in the 1-d case. Think of \code{alpha} as
#       fine-tuning of the penalty.}
# \item{Our method returns an upper-triangular matrix \code{d} when \code{k > 1} so
#       that \code{x \%*\% v == u \%*\% d}. Non-zero
#       elements above the diagonal result from non-orthogonality of the \code{v} matrix,
#       providing a simple interpretation of cumulative information, or explained variance
#       in the PCA case, via the singular value decomposition of \code{d \%*\% t(v)}.}
# }
#
# What if you have no idea for values of the argument \code{n} (the desired sparsity)?
# The reference describes a cross-validation and an ad-hoc approach; neither of which are
# in the package yet. Both are prohibitively computationally expensive for matrices with a huge
# number of columns. A future version of this package will include a revised approach to
# automatically selecting a reasonable sparsity constraint.
#
# Compare with the similar but more general functions \code{SPC} and \code{PMD} in the \code{PMA} package
# by Daniela M. Witten, Robert Tibshirani, Sam Gross, and Balasubramanian Narasimhan.
# The \code{PMD} function can compute low-rank regularized matrix decompositions with sparsity penalties
# on both the \code{u} and \code{v} vectors. The \code{ssvd} function is
# similar to the PMD(*, L1) method invocation of \code{PMD} or alternatively the \code{SPC} function.
# Although less general than \code{PMD}(*),
# the \code{ssvd} function can be faster and more memory efficient for the
# basic sparse PCA problem.
# See \url{https://bwlewis.github.io/irlba/ssvd.html} for more information.
#
# (* Note that the s4vd package by Martin Sill and Sebastian Kaiser, \url{https://cran.r-project.org/package=s4vd},
# includes a fast optimized version of a closely related algorithm by Shen, Huang, and Marron, that penalizes
# both \code{u} and \code{v}.)
#
# @references
# \itemize{
#   \item{Shen, Haipeng, and Jianhua Z. Huang. "Sparse principal component analysis via regularized low rank matrix approximation." Journal of multivariate analysis 99.6 (2008): 1015-1034.}
#   \item{Witten, Tibshirani and Hastie (2009) A penalized matrix decomposition, with applications to sparse principal components and canonical correlation analysis. _Biostatistics_ 10(3): 515-534.}
# }
# @examples
#
# set.seed(1)
# u <- matrix(rnorm(200), ncol=1)
# v <- matrix(c(runif(50, min=0.1), rep(0,250)), ncol=1)
# u <- u / drop(sqrt(crossprod(u)))
# v <- v / drop(sqrt(crossprod(v)))
# x <- u %*% t(v) + 0.001 * matrix(rnorm(200*300), ncol=300)
# s <- ssvd(x, n=50)
# table(actual=v[, 1] != 0, estimated=s$v[, 1] != 0)
# oldpar <- par(mfrow=c(2, 1))
# plot(u, cex=2, main="u (black circles), Estimated u (blue discs)")
# points(s$u, pch=19, col=4)
# plot(v, cex=2, main="v (black circles), Estimated v (blue discs)")
# points(s$v, pch=19, col=4)
#
# # Let's consider a trivial rank-2 example (k=2) with noise. Like the
# # last example, we know the exact number of nonzero elements in each
# # solution vector of the noise-free matrix. Note the application of
# # different sparsity constraints on each column of the estimated v.
# # Also, the decomposition is unique only up to sign, which we adjust
# # for below.
# set.seed(1)
# u <- qr.Q(qr(matrix(rnorm(400), ncol=2)))
# v <- matrix(0, ncol=2, nrow=300)
# v[sample(300, 15), 1] <- runif(15, min=0.1)
# v[sample(300, 50), 2] <- runif(50, min=0.1)
# v <- qr.Q(qr(v))
# x <- u %*% (c(2, 1) * t(v)) + .001 * matrix(rnorm(200 * 300), 200)
# s <- ssvd(x, k=2, n=colSums(v != 0))
#
# # Compare actual and estimated vectors (adjusting for sign):
# s$u <- sign(u) * abs(s$u)
# s$v <- sign(v) * abs(s$v)
# table(actual=v[, 1] != 0, estimated=s$v[, 1] != 0)
# table(actual=v[, 2] != 0, estimated=s$v[, 2] != 0)
# plot(v[, 1], cex=2, main="True v1 (black circles), Estimated v1 (blue discs)")
# points(s$v[, 1], pch=19, col=4)
# plot(v[, 2], cex=2, main="True v2 (black circles), Estimated v2 (blue discs)")
# points(s$v[, 2], pch=19, col=4)
# par(oldpar)
#
ssvd <- function(x, k=1, n=2, maxit=500, tol=1e-3, center=FALSE, scale.=FALSE, alpha=0, tsvd=NULL, ...)
{
  if (alpha < 0  || alpha >= 1) stop("0 <= alpha < 1")
  if (is.logical(center) && center) center <- colMeans(x)
  if (is.logical(scale.))
  {
    if (scale.)
    {
      if (is.numeric(center))
      {
        f <- function(i) sqrt(sum((x[, i] - center[i]) ^ 2) / (nrow(x) - 1L))
        scale. <- vapply(seq(ncol(x)), f, pi, USE.NAMES=FALSE)
      } else scale. <- apply(x, 2L, function(v) sqrt(sum(v ^ 2) / max(1, length(v) - 1L)))
    }
  }
  if (all(n > ncol(x) - 1))
  {
    warning("no sparsity constraints specified")
    return(irlba__irlba(x, k, ...))
  }
  n <- ncol(x) - n
  if (length(n) != k) n <- rep(n, length.out=k) # warn?
  s <- tsvd
  if (is.null(tsvd)) s <- irlba__irlba(x, k, scale=scale., center=center, ...)
  lambda <- c()
  soft <- function(x, u, p)
  {
    y <- crossprod(x, u)
    if (is.numeric(center)) y <- y - sum(u) * center
    if (is.numeric(scale.)) y <- y / scale.
    # apply a column-wise penalty
    a <- abs(y)
    z <- apply(a, 2, sort)
    lambda <<- vapply(seq(length(p)), function(j) (1 - alpha) * z[p[j], j] + alpha * z[p[j] + 1, j], pi, USE.NAMES=FALSE)
    sign(y) * pmax(sweep(a, 2, lambda, `-`), 0)
  }
  s$v <- s$d * s$v
  iter <- 0
  delta_u <- Inf
  while (delta_u > tol && iter < maxit)
  {
    u <- s$u
    s$v <- soft(x, s$u, n)
    if (is.numeric(scale.)) s$v <- s$v / scale.
    if (is.numeric(center))
    {
      xsv <- x %*% s$v - drop(crossprod(center, s$v))
      s$u <- qr.Q(qr(xsv))
    } else
    {
      xsv <- x %*% s$v
      s$u <- qr.Q(qr(xsv))
    }
    # Maintain sign (possibly messed up by QR)
    s$u <- sweep(s$u, 2, apply(xsv, 2, function(x) sign(head(x[x!=0], 1))) / apply(s$u, 2, function(x) sign(head(x[x!=0], 1))), `*`)
    delta_u <- max(1 - diag(abs(crossprod(u, s$u))))
    iter <- iter + 1
  }
  if (iter >= maxit)
    warning("Maximum number of iterations reached before convergence: solution may not be optimal. Consider increasing 'maxit'.")
  s$v <- s$v %*% diag(1 / sqrt(apply(s$v, 2, crossprod)), ncol(s$v), ncol(s$v))
  d <- s$v
  if (is.numeric(scale.)) d <- d / scale.
  d1 <- x %*% d
  if (is.numeric(center)) d1 <- d1 - drop(crossprod(center, d))
  d <- crossprod(s$u, d1)
  list(u = s$u, v = s$v, d = d, iter = iter, lambda = lambda, center=center, scale=scale., n=n, alpha=alpha)
}

# Find a few approximate largest singular values and corresponding
# singular vectors of a matrix.
#
# The randomized method for truncated SVD by P. G. Martinsson and colleagues
# finds a few approximate largest singular values and corresponding
# singular vectors of a sparse or dense matrix. It is a fast and
# memory-efficient way to compute a partial SVD, similar in performance
# for many problems to \code{\link{irlba}}. The \code{svdr} method
# is a block method and may produce more accurate estimations with
# less work for problems with clustered large singular values (see
# the examples). In other problems, \code{irlba} may exhibit faster
# convergence.
#
# Also see an alternate implementation (\code{rsvd}) of this method by N. Benjamin Erichson
# in the https://cran.r-project.org/package=rsvd package.
#
# @param x numeric real- or complex-valued matrix or real-valued sparse matrix.
# @param k dimension of subspace to estimate (number of approximate singular values to compute).
# @param tol stop iteration when the largest absolute relative change in estimated singular
#   values from one iteration to the next falls below this value.
# @param it maximum number of algorithm iterations.
# @param extra number of extra vectors of dimension \code{ncol(x)}, larger values generally improve accuracy (with increased
# computational cost).
# @param center optional column centering vector whose values are implicitly subtracted from each
#   column of \code{A} without explicitly forming the centered matrix (preserving sparsity).
#   Optionally specify \code{center=TRUE} as shorthand for \code{center=colMeans(x)}.
#   Use for efficient principal components computation.
# @param Q optional initial random matrix, defaults to a matrix of size \code{ncol(x)} by \code{k + extra} with
# entries sampled from a normal random distribution.
# @param return.Q if \code{TRUE} return the \code{Q} matrix for restarting (see examples).
# @return
# Returns a list with entries:
# \describe{
#   \item{d:}{ k approximate singular values}
#   \item{u:}{ k approximate left singular vectors}
#   \item{v:}{ k approximate right singular vectors}
#   \item{mprod:}{ total number of matrix products carried out}
#   \item{Q:}{ optional subspace matrix (when \code{return.Q=TRUE})}
# }
# @seealso \code{\link{irlba}}, \code{\link{svd}}, \code{rsvd} in the rsvd package
# @references
# Finding structure with randomness: Stochastic algorithms for constructing
# approximate matrix decompositions N. Halko, P. G. Martinsson, J. Tropp. Sep. 2009.
# @examples
# set.seed(1)
#
# A <- matrix(runif(400), nrow=20)
# svdr(A, 3)$d
#
# # Compare with svd
# svd(A)$d[1:3]
#
# # Compare with irlba
# irlba(A, 3)$d
#
# \dontrun{
# # A problem with clustered large singular values where svdr out-performs irlba.
# tprolate <- function(n, w=0.25)
# {
#   a <- rep(0, n)
#   a[1] <- 2 * w
#   a[2:n] <- sin( 2 * pi * w * (1:(n-1)) ) / ( pi * (1:(n-1)) )
#   toeplitz(a)
# }
#
# x <- tprolate(512)
# set.seed(1)
# tL <- system.time(L <- irlba(x, 20))
# tR <- system.time(R <- svdr(x, 20))
# S <- svd(x)
# plot(S$d)
# data.frame(time=c(tL[3], tR[3]),
#            error=sqrt(c(crossprod(L$d - S$d[1:20]), crossprod(R$d - S$d[1:20]))),
#            row.names=c("IRLBA", "Randomized SVD"))
#
# # But, here is a similar problem with clustered singular values where svdr
# # doesn't out-perform irlba as easily...clusters of singular values are,
# # in general, very hard to deal with!
# # (This example based on https://github.com/bwlewis/irlba/issues/16.)
# set.seed(1)
# s <- svd(matrix(rnorm(200 * 200), 200))
# x <- s$u %*% (c(exp(-(1:100)^0.3) * 1e-12 + 1, rep(0.5, 100)) * t(s$v))
# tL <- system.time(L <- irlba(x, 5))
# tR <- system.time(R <- svdr(x, 5))
# S <- svd(x)
# plot(S$d)
# data.frame(time=c(tL[3], tR[3]),
#            error=sqrt(c(crossprod(L$d - S$d[1:5]), crossprod(R$d - S$d[1:5]))),
#            row.names=c("IRLBA", "Randomized SVD"))
# }
irlba__svdr <- function(x, k, tol=1e-5, it=100L, extra=min(10L, dim(x) - k), center=NULL, Q=NULL, return.Q=FALSE)
{
  eps2 <- .Machine$double.eps ^ (4 / 5)
  n <- min(ncol(x), k + extra)
  if (isTRUE(center)) center <- colMeans(x)
  if (is.null(Q)) Q <-  matrix(stats::rnorm(ncol(x) * n), ncol(x))
  d <- rep(0, k)
  for (j in 1:it)
  {
    if (is.null(center))
    {
      Q <- qr.Q(qr(x %*% Q))
      B <- crossprod(x, Q)
      Q <- qr.Q(qr(B))
    } else
    {
      Q <- qr.Q(qr(x %*% Q - rep(1, nrow(x)) %*% crossprod(center, Q)))
      B <- crossprod(Q, x) - tcrossprod(crossprod(Q, rep(1, nrow(x))), center)
      Q <- qr.Q(qr(t(B)))
    }
    d1 <- svd(B, nu=0, nv=0)$d[1:k]
    idx <- d1 > eps2
    if (all(! idx)) break
    if (max(abs((d1[idx] - d[idx]) / d[idx])) < tol) break
    d <- d1
  }
  if (return.Q) Q1 <- Q
  if (is.null(center))
  {
    Q <- qr.Q(qr(x %*% Q))
    B <- crossprod(Q, x)
  } else
  {
    Q <- qr.Q(qr(x %*% Q - rep(1, nrow(x)) %*% crossprod(center, Q)))
    B <- crossprod(Q, x) - tcrossprod(crossprod(Q, rep(1, nrow(x))), center)
  }
  s <- svd(B)
  s$u <- Q %*% s$u
  s$u <- s$u[, 1:k]
  s$d <- s$d[1:k]
  s$v <- s$v[, 1:k]
  s$mprod <- 2 * j + 1
  if (return.Q) s$Q <- Q1
  s
}

# ---------------------------------------------------------------------
# Internal supporting functions
# ---------------------------------------------------------------------
# General real/complex crossprod
cross <- function(x, y)
{
  if (missing(y))
  {
    if (is.complex(x)) return(abs(Conj(t(x)) %*% x))
    return(crossprod(x))
  }
  if (!is.complex(x) && !is.complex(y)) return(crossprod(x, y))
  Conj(t(x)) %*% y
}

# Euclidean norm
norm2 <- function(x)
{
  drop(sqrt(cross(x)))
}

# Orthogonalize vectors Y against vectors X.
orthog <- function(Y, X)
{
  dx2 <- dim(X)[2]
  if (is.null(dx2)) dx2 <- 1
  dy2 <- dim(Y)[2]
  if (is.null(dy2)) dy2 <- 1
  if (dx2 < dy2) doty <- cross(X, Y)
  else doty <- Conj(t(cross(Y, X)))
  Y - X %*% doty
}

# Convergence tests
# Input parameters
# Bsz            Number of rows of the bidiagonal matrix B
# tol
# k_org
# Bsvd           svd list of small matrix B
# residuals
# k
# Smax
# lastsv, svtol, work, S
#
# Output parameter list
# converged      TRUE/FALSE
# k              Number of singular vectors returned
convtests <- function(Bsz, tol, k_org, Bsvd, residuals, k, Smax, lastsv, svtol, maxritz, work, S)
{
  # Converged singular triplets
  subspace_converged <- residuals[1:k_org] < tol * Smax
  # Converged fixed point triplets
  if (is.null(lastsv)) lastsv <- 0 * Bsvd$d
  delta_converged <- (abs(Bsvd$d[1:k_org] - lastsv[1:k_org]) / Bsvd$d[1:k_org])  < svtol
  len_res <- sum(subspace_converged & delta_converged) # both
  if (is.na(len_res)) len_res <- 0
  if (len_res >= k_org) return(list(converged=TRUE, k=k))
  if (S == 0) return(list(converged=TRUE, k=k))
  # Not converged yet...
  # Adjust k to include more vectors as the number of vectors converge, but not
  # too many (maxritz):
  augment <- min(sum(subspace_converged), maxritz)
  k <- min(max(k, k_org + augment), work - 1)
  list(converged=FALSE, k=k)
}

message_once <- function(..., flag)
{
  if (flag$flag) return()
  flag$flag <- TRUE
  message(...)
}
