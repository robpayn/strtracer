#' calculate the uptakes for provided concentrations based on a saturating hyperbolic kinetic model
#' 
#' @export
hyperbolic <- function(umax, halfsat, conc)
{
  return((umax * conc) / (halfsat + conc));
}

#' calculate the net uptakes for provided concentrations based on a saturating hyperbolic kinetic model
#' 
#' @export
hyperbolicnet <- function(umax, halfsat, concadd, concbkg)
{
  return(hyperbolic(umax, halfsat, (concadd + concbkg)) - hyperbolic(umax, halfsat, concbkg));
}

#' Create a graphics device based on a provided device name
#' 
#' @export
createDevice <- function(deviceName, width, height) 
{
   switch(
      deviceName, 
      windows = windows(width = width, height = height),
      message("Using active graphics device.")
      );
}

#' Create a blank plot window
#' 
#' @export
createBlankPlot <- function(
   x = 1,
   y = 1,
   type = "n",
   xlim, 
   ylim, 
   xlab, 
   ylab, 
   ...
   )
{
   plot(
      x = x,
      y = y,
      xlim = xlim,
      ylim = ylim,
      type = type,
      xlab = xlab,
      ylab = ylab,
      ...
      );
}
