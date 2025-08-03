stime <- function() {
  format(Sys.time(), "%T")
}

# message with a time stamp
# appears only if called from an environment where a logical verbose = TRUE
# OR force = TRUE
tsmessage <- function(..., domain = NULL, appendLF = TRUE, force = FALSE,
                      time_stamp = TRUE) {
  verbose <- get0("verbose", envir = sys.parent())

  if (force || (!is.null(verbose) && verbose)) {
    msg <- ""
    if (time_stamp) {
      msg <- paste0(stime(), " ")
    }
    message(msg, ..., domain = domain, appendLF = appendLF)
    utils::flush.console()
  }
}

# log vector information
summarize <- function(X, msg = "") {
  summary_X <- summary(X, digits = max(3, getOption("digits") - 3))
  tsmessage(msg, ": ", paste(names(summary_X), ":", summary_X, "|",
    collapse = ""
  ),
  force = get0("verbose", envir = sys.parent())
  )
}

# pluralize("thread", 1) => "1 thread"
# pluralize("thread", 2) => "2 threads"
pluralize <- function(str, n, prefix = NULL, inc_num = TRUE) {
  if (n == 0) {
    return("")
  }
  ret <- paste0(str, ifelse(n != 1, "s", ""))
  if (inc_num) {
    ret <- paste0(n, " ", ret)
  }
  if (!is.null(prefix)) {
    ret <- paste0(prefix, " ", ret)
  }
  ret
}

# convert data frame to matrix using numeric columns
x2m <- function(X) {
  if (!methods::is(X, "matrix")) {
    m <- as.matrix(X[, which(vapply(X, is.numeric, logical(1)))])
  } else {
    m <- X
  }
  m
}

# given a metric argument, returns a list containing:
# metrics - the input list with any members called "categorical" removed
# categoricals - a vector of the categorical ids
find_categoricals <- function(metrics) {
  res <- list(
    metrics = metrics
  )
  if (is.list(metrics)) {
    cat_pos <- grep("categorical", names(metrics))
    if (length(cat_pos) > 0) {
      cat_ids <- unlist(metrics[cat_pos])
      names(cat_ids) <- NULL
      res <- list(
        metrics = metrics[-cat_pos],
        categoricals = cat_ids
      )
    }
  }
  res
}

# Splits a list into its named and unnamed components:
# > lsplit_unnamed(list(1:10, pca_center = FALSE))
# $named
# $named$pca_center
# [1] FALSE
#
#
# $unnamed
# $unnamed[[1]]
#  [1]  1  2  3  4  5  6  7  8  9 10
lsplit_unnamed <- function(l) {
  lnames <- names(l)
  if (is.null(lnames)) {
    return(list(unnamed = l))
  }
  is_named <- lnames != ""
  nids <- which(is_named)
  uids <- which(!is_named)
  if (length(uids) == 0) {
    return(list(named = l[nids]))
  }

  list(
    named = l[nids],
    unnamed = l[uids]
  )
}

# Do work and update a progress bar
progress_for <- function(n, nchunks, fun) {
  message("0%   10   20   30   40   50   60   70   80   90   100%")
  message("[----|----|----|----|----|----|----|----|----|----|")
  remaining <- n
  chunk_end <- 0
  for (i in 1:nchunks) {
    chunk_start <- chunk_end + 1
    chunk_end <- chunk_start + round(remaining / (nchunks - i + 1)) - 1
    remaining <- remaining - (chunk_end - chunk_start + 1)
    fun(chunk_start, chunk_end)

    message("*", appendLF = FALSE)
    utils::flush.console()
  }
  message("|")
}

checkna <- function(X) {
  if (!is.null(X) && any(is.na(X))) {
    stop("Missing values found in 'X'")
  }
}

check_graph <- function(graph, expected_rows = NULL, expected_cols = NULL,
                        bipartite = FALSE) {
  idx <- graph$idx
  dist <- graph$dist
  if (!methods::is(idx, "matrix")) {
    stop("neighbor graph must contain an 'idx' matrix")
  }
  if (!methods::is(dist, "matrix")) {
    stop("neighbor graph must contain a 'dist' matrix")
  }
  if (!all(dim(idx) == dim(dist))) {
    stop("'idx' and 'dist' matrices must have identical dimensions")
  }
  # graph may be our only source of input data, in which case no other source
  # to validate from
  if (!is.null(expected_rows)) {
    if (nrow(idx) != expected_rows) {
      stop("idx matrix has unexpected number of rows")
    }
  }
  if (!is.null(expected_cols) && !is.na(expected_cols)) {
    if (ncol(idx) != expected_cols) {
      stop("idx matrix has unexpected number of columns")
    }
  }
  # graph should not contain missing data
  if (any(is.na(idx)) || any(is.na(dist))) {
    stop("neighbor graph contains missing data")
  }

  # if looking at neighbors within one graph there can't be more neighbors
  # than observations
  if (!bipartite) {
    if (ncol(idx) > nrow(idx)) {
      stop("Invalid neighbors: number exceeds number of observations")
    }
    if (max(idx) > nrow(idx)) {
      stop("Invalid neighbors: max index exceeds number of observations")
    }
  }
}

check_sparse_graph <- function(graph, expected_rows = NULL,
                               expected_cols = NULL, bipartite = FALSE) {
  if (!is.null(expected_rows)) {
    if (nrow(graph) != expected_rows) {
      stop("Sparse distance matrix has unexpected number of rows")
    }
  }
  if (!is.null(expected_cols)) {
    if (ncol(graph) != expected_cols) {
      stop("Sparse distance matrix has unexpected number of cols")
    }
  }
  if (!bipartite) {
    if (nrow(graph) != ncol(graph)) {
      stop("Sparse distance matrix must have same number of rows and cols")
    }
  }
}

check_graph_list <- function(graph_list, expected_rows = NULL,
                             expected_cols = NULL, bipartite = FALSE) {
  if (nn_is_single(graph_list)) {
    graph_list <- list(graph_list)
  }
  num_nns <- length(graph_list)
  if (num_nns == 0) {
    stop("precalculated graph list is empty")
  }
  for (i in 1:num_nns) {
    graph <- graph_list[[i]]
    if (is.list(graph)) {
      check_graph(graph, expected_rows, expected_cols, bipartite = bipartite)
    } else if (is_sparse_matrix(graph)) {
      check_sparse_graph(graph, expected_rows, expected_cols,
        bipartite = bipartite
      )
    } else {
      stop("Unknown neighbor data format")
    }
  }
  num_nns
}

nn_graph_row_names_list <- function(graph_list) {
  if (nn_is_single(graph_list)) {
    graph_list <- list(graph_list)
  }
  xnames <- NULL
  for (i in 1:length(graph_list)) {
    graph <- graph_list[[i]]
    if (is.list(graph)) {
      xnames <- nn_graph_row_names(graph)
    } else if (is_sparse_matrix(graph)) {
      xnames <- row.names(graph)
    } else {
      stop("Unknown neighbor data format")
    }
    if (!is.null(xnames)) {
      break
    }
  }
  xnames
}

# from a nn graph (or list) get the first non-NULL row names
nn_graph_row_names <- function(graph) {
  xnames <- NULL
  if (!is.null(row.names(graph$idx))) {
    xnames <- row.names(graph$idx)
  }
  if (is.null(xnames) && !is.null(row.names(graph$dist))) {
    xnames <- row.names(graph$dist)
  }
  xnames
}

nn_graph_nbrs_list <- function(graph_list) {
  if (nn_is_single(graph_list)) {
    graph_list <- list(graph_list)
  }
  sapply(graph_list, nn_graph_nbrs)
}

# from a nn graph (or list) get the number of neighbors
nn_graph_nbrs <- function(graph) {
  if (is.list(graph)) {
    ncol(graph$idx)
  } else if (is_sparse_matrix(graph)) {
    NA
  } else {
    stop("Unknown neighbor data format")
  }
}

is_sparse_matrix <- function(m) {
  methods::is(m, "sparseMatrix")
}

# Add the (named) values in l2 to l1.
# Use to override default values in l1 with user-supplied values in l2
lmerge <- function(l1, l2) {
  for (name in names(l2)) {
    l1[[name]] <- l2[[name]]
  }
  l1
}

range_scale <- function(x, min = 0, max = 1) {
  (x - min(x)) / (max(x) - min(x)) * (max - min) + min
}

is_installed <- function(pkgname) {
  requireNamespace(pkgname,
    quietly = TRUE,
  )
  isNamespaceLoaded(pkgname)
}

is_win7 <- function() {
  sys_info <- Sys.info()
  sys_info[["sysname"]] == "Windows" &&
    strsplit(sys_info["release"], split = " ")$release[[1]] == "7"
}
