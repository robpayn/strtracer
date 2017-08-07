require(XML);

# CLASS Experiment ####

#' Create a new stream tracer experiment object
#' 
#' @export
#' @param discharge Stream discharge during experiment
#' @param streamLength Length of experimental stream reach
#' @param streamWidth Average width of stream channel
#' @param streamDepth Average depth of stream channel
#' @param xSectionArea Cross-sectional area of stream channel.  
#'    Calculated as product of width and depth by default.
#' @param streamVel Depth-averaged velocity of stream channel.  
#'    Calculated as the discharge divided by the cross-sectional area by default.
#' @param travelTime Average transport time of water along the reach.  
#'    Calculated as reach length devided by the depth averaged velocity by default.
#' @param conserveSolute Dataframe with the time series of conservative solute breakthrough concentrations
#' @param activeSolute Dataframe with the time series of active solute breakthrough concentrations
#' @param conserveBkg Background concentration for conservative solute
#' @param activeBkg Background concentration for active solute
#' @param injectRatio Ratio of the active to conservative tracer in injectate.
#'    Default value is unity.
#' @return a new experiment object configured by arguments
Experiment <- function(
   discharge,
   streamLength,
   streamWidth,
   streamDepth,
   xSectionArea = streamWidth * streamDepth,
   streamVel = discharge / (streamWidth * streamDepth),
   travelTime = streamLength / 
      (discharge / (streamWidth * streamDepth)),
   conserveSolute,
   activeSolute,
   conserveBkg,
   activeBkg,
   injectRatio = 1
   )
{
   experiment <- new.env();
   
   experiment$discharge <- discharge;
   experiment$streamLength <- streamLength;
   experiment$streamWidth <- streamWidth;
   experiment$streamDepth <- streamDepth;
   experiment$xSectionArea <- xSectionArea;
   experiment$streamVel <- streamVel;
   experiment$travelTime <- travelTime;
   experiment$conserveSolute <- conserveSolute;
   experiment$activeSolute <- activeSolute;
   experiment$conserveBkg <- conserveBkg;
   experiment$activeBkg <- activeBkg;
   experiment$injectRatio <- injectRatio;
   
   class(experiment) <- c("Experiment", class(experiment));
   return(experiment);
}

#' Plot the conservative tracer breakthrough curve
#' 
#' @export
#' @param experiment Experiment for which the plot is generated
#' @param ... Other parameters
#' @return Output written to graphics device, nothing is returned
plotConservative <- function(experiment, ...)
{
   UseMethod("plotConservative", experiment);
}

#' Plot the conservative tracer breakthrough curve
#' 
#' @export
#' @param experiment Experiment for which the plot is generated
#' @param ... Other parameters
#' @return Output written to graphics device, nothing is returned
plotConservative.Experiment <- function(
   experiment,
   device = "default",
   width = 8,
   height = 6,
   columns = 3:length(experiment$conserveSolute),
   xfactor = 1,
   yfactor = 1,
   xlim = c(
      min(experiment$conserveSolute$Time),
      max(experiment$conserveSolute$Time)
      ),
   ylim = c(
      0,
      max(
         if (backgroundCorrect) (experiment$conserveSolute[,columns] - experiment$conserveBkg)
         else experiment$conserveSolute[,columns]
         )
      ),
   xlab = "Time",
   ylab = "Concentration",
   ratio = 1,
   backgroundCorrect = TRUE,
   ...
   ) 
{
   createDevice(device, width, height);
   par(...);
   createBlankPlot(
      xlim = xlim * xfactor, 
      ylim = ylim * yfactor * ratio, 
      xlab = xlab, 
      ylab = ylab
      );
   for (column in columns)
   {
      lines(
         x = experiment$conserveSolute$Time * xfactor,
         y = if (backgroundCorrect) 
               ((experiment$conserveSolute[[column]] - experiment$conserveBkg) * yfactor * ratio)
            else (experiment$conserveSolute[[column]] * yfactor * ratio)
      )
   }
}

#' Plot the active tracer breakthrough curve
#' 
#' @export
plotActive <- function(experiment, ...)
{
   UseMethod("plotActive", experiment);
}

#' Plot the active tracer breakthrough curve for a tracer experiment
#' 
#' @export
plotActive.Experiment <- function(
   experiment,
   columns = 3:length(experiment$conserveSolute),
   xfactor = 1,
   yfactor = 1,
   ratio = experiment$injectRatio,
   activeColor = "red",
   window = NULL,
   ...
   ) 
{
   plotConservative.Experiment(
      experiment = experiment,
      columns = columns,
      xfactor = xfactor,
      yfactor = yfactor,
      ratio = ratio,
      ...
      );
   for (column in columns)
   {
      lines(
         x = experiment$activeSolute$Time * xfactor,
         y = (experiment$activeSolute[[column]] - experiment$activeBkg) 
            * yfactor,
         col = activeColor
      )
   }
   if (!is.null(window))
   {
      abline(v = window * xfactor, lty = "dashed", col = "red");
   }
}

# CLASS ExperimentSlug ####

#' Create a new instance of an instantaneous release (slug) tracer experiment
#' 
#' @export
ExperimentSlug <- function(
   experiment = Experiment(injectRatio = injectRatio, ...),
   releaseTime,
   conserveMass,
   activeMass,
   injectRatio = activeMass / conserveMass,
   ...
   )
{
   experiment$releaseTime <- releaseTime;
   experiment$conserveMass <- conserveMass;
   experiment$activeMass <- activeMass;
   experiment$injectRatio <- injectRatio;

   class(experiment) <- c("ExperimentSlug", class(experiment));
   return(experiment);   
}

#' Plot the conservative tracer breakthrough curve for a slug experiment
#' 
#' @export
plotConservative.ExperimentSlug <- function(
   experiment,
   xfactor = 1,
   releaseTimeCol = "black",
   releaseTimeLty = "dashed",
   ...
   ) 
{
   plotConservative.Experiment(
      experiment = experiment, 
      xfactor = xfactor, 
      ...
      );
   abline(
      v = experiment$releaseTime * xfactor, 
      lty = releaseTimeLty, 
      col = releaseTimeCol
   );
}

#' Plot the active tracer breakthrough curve for a slug experiment
#' 
#' @export
plotActive.ExperimentSlug <- function(
   experiment,
   xfactor = 1,
   ratio = experiment$injectRatio,
   releaseTimeCol = "black",
   releaseTimeLty = "dashed",
   ...
   ) 
{
   plotActive.Experiment(
      experiment = experiment,
      xfactor = xfactor,
      ratio = ratio,
      ...
      );
   abline(
      v = experiment$releaseTime * xfactor, 
      lty = releaseTimeLty, 
      col = releaseTimeCol
      );
}
