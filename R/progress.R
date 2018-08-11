# An RC class to produce a progress bar like RcppProgress does for C++ code.
# create:
#   progress <- Progress$new(max, display)
# max is the final value (e.g. last index in a loop)
# display is a logical for whether to actually display the bar
# use:
#   progress$increment()
# increments the counter towards max. A new star will be displayed on the bar if sufficient
# progress was made. Nothing happens if the bar was created with display = FALSE.
#
# Typical usage:
# do_stuff_ntimes <- function(n, verbose = FALSE, ...)
#   progress <- Progress$new(max = n, display = verbose)
#   for (i in 1:n) {
#     do_more_stuff(i, ...)
#     progress$increment()
#   }
# }
#
Progress <- setRefClass("Progress",
  fields = list(
    value = "numeric",
    max = "numeric",
    curr_stars = "numeric",
    max_stars = "numeric",
    display = "logical"
  ),
  methods = list(
    initialize = function(max, display = TRUE) {
      max_stars <<- 51 # length of the progress bar
      value <<- 0
      curr_stars <<- 0
      max <<- max
      display <<- display
      if (display) {
        message("0%   10   20   30   40   50   60   70   80   90   100%")
        message("[----|----|----|----|----|----|----|----|----|----|")
      }
    },
    increment = function() {
      if (display && curr_stars < max_stars) {
        value <<- value + 1
        num_stars <- round(max_stars * value / max)
        if (num_stars > curr_stars) {
          # Number of new stars to print
          num_new_stars <- num_stars - curr_stars

          # If we are going to reach the end of the progress bar
          # save space for the terminal "|"
          if (curr_stars + num_new_stars >= max_stars) {
            num_new_stars <- num_new_stars - 1
          }
          new_stars <- paste(rep("*", num_new_stars), collapse = "")

          message(new_stars, appendLF = FALSE)
          flush.console()
          curr_stars <<- num_stars
        }
        if (curr_stars >= max_stars) {
          # The terminal "|" character that appears instead of a *
          message("|")
        }
      }
    }
  )
)
