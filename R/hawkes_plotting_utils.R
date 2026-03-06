#' Create a Modern Base R Canvas
#' @description Sets up a clean plot area with modern margins, tight axes, and a
#'   soft grid.
#' @param xlim Numeric vector of length 2 (e.g., c(0, T)).
#' @param ylim Numeric vector of length 2 (e.g., c(0, max_y)).
#' @param title String for the plot title.
#' @param xlab String for the x-axis label. Default "Time (t)".
#' @param ylab String for the y-axis label. Default intensity symbol.
#' 
#' @importFrom graphics par grid
#' @export
plot_base_canvas <- function(
  xlim, 
  ylim, 
  title, 
  xlab = "Time (t)", 
  ylab = expression(lambda^"*"*(t))
) {
  # Modern styling: 
    # bty="l" (L-shaped), 
    # las=1 (horizontal text), 
    # yaxs="i" (tight axes)
  par(mar = c(4, 4, 3, 1), bty = "l", las = 1, family = "sans")
  
  plot(0, type = "n", xlim = xlim, ylim = ylim, 
       xaxs = "i", yaxs = "i",
       xlab = xlab, ylab = ylab, 
       main = title, col.main = "#2E3440", cex.main = 1.1)
  
  # Soft horizontal grid lines
  grid(nx = NA, ny = NULL, col = "gray93", lty = "solid")
}

#' Add Intensity Line to Plot
#' @description Layers a muted blue intensity curve onto an existing plot.
#' @param t Numeric vector of time points (the grid).
#' @param intensities Numeric vector of intensity values.
#' @param col Hex code for the line colour. Default is Nord Frost Blue.
#' 
#' @importFrom graphics lines
#' @export
add_intensity <- function(t, intensities, col = "#5E81AC") {
  lines(t, intensities, col = col, lwd = 1.8)
}

#' Add Event Ticks to Plot
#' @description Adds flat-ended rug ticks representing event times.
#' @param H_t Numeric vector of event timestamps.
#' @param col Hex code for the tick colour. Default is Nord Aurora Red.
#' 
#' @importFrom graphics rug
#' @export
add_events <- function(H_t, col = "#BF616A") {
  # lend=1 ensures flat, professional-looking bar ends
  rug(H_t, col = col, lwd = 1, ticksize = 0.03, lend = 1)
}

#' Add Counting Process (Step Plot)
#' @description Layers a step-function representing the counting process N(t).
#' @param H_t Numeric vector of event timestamps.
#' @param T_max The end of the observation window.
#' @param col Hex code for the step line. Default is Nord Slate.
#' 
#' @importFrom graphics lines
#' @export
add_counting_process <- function(H_t, T_max, col = "#4C566A") {
  # We start at (0,0), then jump at each event time
  # x-coords: 0 followed by the event times
  # y-coords: 0 followed by 1, 2, ..., n
  x_vals <- c(0, H_t)
  y_vals <- 0:length(H_t)
  
  # type = "s" creates the step starting from the left
  lines(x_vals, y_vals, type = "s", col = col, lwd = 1.8)
}
