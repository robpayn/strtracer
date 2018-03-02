setOldClass(Classes = c("Simulation", "Experiment"));

# Class HyperbolicAnalysis definition and constructor ####

#' Class HyperbolicAnalysis
#' 
#' Class \code{HyperbolicAnalysis} defines an analysis of whole-stream
#' kinetic behavior based on solute uptake experiments
#' 
#' @export HyperbolicAnalysis
#' @exportClass HyperbolicAnalysis

HyperbolicAnalysis <- setClass(
   Class = "HyperbolicAnalysis",
   contains = "environment",
   );

#' Constructor method for HyperbolicAnalysis class.
#' 
#' @rdname HyperbolicAnalysis-class
#' @usage S4 Constructors:
#' HyperbolicAnalysis(experiment, time, indices)
#' new(Class, experiment, time, indices)
#' @param Class Should be designated "HyperbolicAnalysis" to 
#'      instantiate this class
#' @param experiment Tracer experiment object
#' @param time Vector of times associated with solute 
#'      uptake estimates
#' @param indices Vector of indices used to extract
#'      appropriate solute concentrations
#' @return The object of class \code{HyperbolicAnalysis} created
#'      by the constructor

setMethod(
   f = "initialize",
   signature = "HyperbolicAnalysis",
   definition = function(
      .Object, 
      experiment,
      time,
      indices
      )
      {
         analysis <- callNextMethod(.Object);
         
         analysis$experiment <- experiment;
         analysis$indices <- indices;
         
         metricsLength <- length(time);
         analysis$metrics <- data.frame(
            time = time,
            conserve = numeric(metricsLength),
            active = numeric(metricsLength),
            conservebc = numeric(metricsLength),
            activebc = numeric(metricsLength),
            activenr = numeric(metricsLength),
            travelTime = numeric(metricsLength),
            ceffinject = numeric(metricsLength),
            cefftot = numeric(metricsLength),
            k = numeric(metricsLength),
            sw = numeric(metricsLength),
            vf = numeric(metricsLength),
            u = numeric(metricsLength),
            stringsAsFactors = FALSE
            );
         
         return(analysis);
      }
   );

# HyperbolicAnalysis.calcMetrics method ####

#' Generic method calcMetrics
#' 
#' @exportMethod calcMetrics

setGeneric(
   name = "calcMetrics",
   def = function(analysis, ...) { standardGeneric("calcMetrics") }
   );

#' Calculate the solute uptake metrics for a 
#' hyperbolic saturating kinetic analysis
#'
#' @aliases calcMetrics.HyperbolicAnalysis
#' @exportMethod calcMetrics

setMethod(
   f = "calcMetrics",
   signature = "HyperbolicAnalysis",
   definition = function(analysis, ...)
      {
      analysis$metrics$cefftot <- 
         analysis$experiment$activeBkg + analysis$metrics$ceffinject;
      analysis$metrics$sw <- 
         analysis$experiment$streamVel / analysis$metrics$k;
      analysis$metrics$vf <- 
         analysis$metrics$k * analysis$experiment$streamDepth;
      analysis$metrics$u <- 
         analysis$metrics$vf * analysis$metrics$ceffinject;
      }
   );

# HyperbolicAnalysis.runAll method ####

#' Generic method runAll
#' 
#' @exportMethod runAll

setGeneric(
   name = "runAll",
   def = function(analysis, ...) { standardGeneric("runAll") }
   );

#' Run all forms of the hyperbolic saturating kinetic analysis
#' 
#' @aliases runAll.HyperbolicAnalysis
#' @exportMethod runAll

setMethod(
   f = "runAll",
   signature = "HyperbolicAnalysis",
   definition = function(analysis, ...) 
      {
         runSw(analysis, ...);
         runU(analysis, ...);
         runVf(analysis, ...);
         runTau(analysis, ...);
      }
   );

# HyperbolicAnalysis.runSw method ####

#' Method runSw
#' 
#' @exportMethod runSw

setGeneric(
   name = "runSw",
   def = function(analysis) { standardGeneric("runSw") }
   );

#' Run the uptake length hyperbolic saturating kinetic analysis
#'
#' @aliases runSw.HyperbolicAnalysis
#' @exportMethod runSw

setMethod(
   f = "runSw",
   signature = "HyperbolicAnalysis",
   definition = function(analysis) 
      {
         # Regression of added solute uptake length vs. concentration
         lmresults <- lm(
           sw ~ cefftot,
           data = analysis$metrics
         );
         analysis$swEstimates <- list(
            intercept = as.numeric(lmresults$coefficients["(Intercept)"]),
            slope = as.numeric(lmresults$coefficients["cefftot"]),
            swamb = numeric(length = 1),
            umax = numeric(length = 1),
            halfsat = numeric(length = 1),
            uamb = numeric(length = 1)
            );
         analysis$swEstimates$swamb <- analysis$swEstimates$intercept;
         analysis$swEstimates$halfsat <- 
            ( (analysis$swEstimates$intercept + analysis$swEstimates$slope * analysis$experiment$activeBkg) 
            / analysis$swEstimates$slope ) - analysis$experiment$activeBkg;
         analysis$swEstimates$umax <- ( analysis$experiment$discharge * 
            (analysis$swEstimates$halfsat + analysis$experiment$activeBkg) ) /
            (analysis$experiment$streamWidth * analysis$swEstimates$slope * analysis$swEstimates$halfsat)
         analysis$swEstimates$uamb <- (analysis$experiment$discharge * analysis$experiment$activeBkg) /
               (analysis$experiment$streamWidth * analysis$swEstimates$swamb);
      }
   );

# HyperbolicAnalysis.runU method ####

#' Method runU
#' 
#' @exportMethod runU

setGeneric(
   name = "runU",
   def = function(analysis, ...) { standardGeneric("runU") }
   );

#' Run the uptake nonlinear hyperbolic saturating kinetic analysis
#'
#' @aliases runU.HyperbolicAnalysis
#' @exportMethod runU

setMethod(
   f = "runU",
   signature = "HyperbolicAnalysis",
   definition = function(
      analysis,
      fixedParameters = list(umaxp = NULL, halfsatp = NULL),
      initialEstimates = NULL
      ) 
      {
         # Nonlinear regression of hyperbolic function for uptake vs. concentration
         if (is.null(initialEstimates))
         {
            start = list(
               umaxp = analysis$swEstimates$umax,
               halfsatp = analysis$swEstimates$halfsat
               );
         }
         else
         {
            start = initialEstimates;
         }
         activeBkg <- analysis$experiment$activeBkg;
         
         if (!is.null(fixedParameters$umaxp)) umaxp <- fixedParameters$umaxp;
         if (!is.null(fixedParameters$halfsatp)) halfsatp <- fixedParameters$halfsatp;
         nlsresults <- nls(
            u ~ hyperbolicnet(
              umax = umaxp, 
              halfsat = halfsatp, 
              concadd = ceffinject, 
              concbkg = activeBkg
              ),
            data = analysis$metrics,
            start = start
            );
         if (is.null(fixedParameters$umaxp))
         {
            umaxest <- summary(nlsresults)$coefficients["umaxp","Estimate"];
         }
         else
         {
            umaxest <- fixedParameters$umaxp;
         }
         if (is.null(fixedParameters$halfsatp))
         {
            halfsatest <- summary(nlsresults)$coefficients["halfsatp","Estimate"];
         }
         else
         {
            halfsatest <- fixedParameters$halfsat;
         }
         analysis$uEstimates <- list(
            umax = umaxest,
            halfsat = halfsatest,
            uamb = hyperbolic(
               umax = umaxest, 
               halfsat = halfsatest, 
               conc = analysis$experiment$activeBkg
               )
            );
      }
   );

# HyperbolicAnalysis.runVf method ####

#' Method runVf
#' 
#' @exportMethod runVf

setGeneric(
   name = "runVf",
   def = function(analysis) { standardGeneric("runVf") }
   );

#' Run the uptake velocity hyperbolic saturating kinetic analysis
#'
#' @aliases runVf.HyperbolicAnalysis
#' @exportMethod runVf

setMethod(
   f = "runVf",
   signature = "HyperbolicAnalysis",
   definition = function(analysis) 
      {
         ineff <- 1 / analysis$metrics$vf;
         lmresults <- lm(
            ineff ~ cefftot,
            data = analysis$metrics
         );
         intercept = as.numeric(lmresults$coefficients["(Intercept)"]);
         slope = as.numeric(lmresults$coefficients["cefftot"]);
         halfsat = intercept / slope;
         vfambest = 1 / intercept;
         analysis$vfEstimates = list(
            intercept = intercept,
            slope = slope,
            vfamb = vfambest,
            umax = (halfsat + analysis$experiment$activeBkg) /
               intercept,
            halfsat = halfsat,
            uamb = vfambest * analysis$experiment$activeBkg
         );
      }
   );

# HyperbolicAnalysis.runTau method ####

#' Method runTau
#' 
#' @exportMethod runTau

setGeneric(
   name = "runTau",
   def = function(analysis) { standardGeneric("runTau") }
   );

#' Run the reaction time hyperbolic saturating kinetic analysis
#'
#' @aliases runTau.HyperbolicAnalysis
#' @exportMethod runTau

setMethod(
   f = "runTau",
   signature = "HyperbolicAnalysis",
   definition = function(analysis) 
      {
         tau <- 1 / analysis$metrics$k;
         lmresults <- lm(
            tau ~ cefftot,
            data = analysis$metrics
         );
         tauAmb = as.numeric(lmresults$coefficients["(Intercept)"]);
         slope = as.numeric(lmresults$coefficients["cefftot"]);
         halfsat = tauAmb / slope;
         analysis$tauEstimates = list(
            intercept = tauAmb,
            slope = slope,
            tauAmb = tauAmb,
            umax = analysis$experiment$streamDepth * 
               (halfsat + analysis$experiment$activeBkg) / tauAmb,
            halfsat = halfsat,
            uamb = (analysis$experiment$streamDepth * analysis$experiment$activeBkg) / 
               tauAmb
         );
      }
   );

# HyperbolicAnalysis.plot method ####

#' Method plot
#' 
#' @exportMethod plot

setGeneric(name = "plot");

#' Plot a basic visualization of the analyses
#'
#' @aliases plot.HyperbolicAnalysis
#' @exportMethod plot

setMethod(
   f = "plot",
   signature = "HyperbolicAnalysis",
   definition = function(
      x, 
      y = NA, 
      xlab = "Effective Concentration", 
      xlab.cex = 1,
      vf.ylab = bquote(paste(v[f])),
      u.ylab = "U",
      sw.ylab = bquote(paste(s[w])),
      tau.ylab = bquote(paste(tau[NET])),
      actual = TRUE,
      actual.col = "red",
      actual.lty = "dashed",
      ...
      ) 
      {
         par(
            mar = c(2, 4.5, 1, 1),
            oma = c(2, 0, 0, 0),
            mfrow = c(2, 2),
            ...
            );
         plotVf(analysis = x, xlab = "", ylab = vf.ylab);
         if (actual) {
            linesVfModel(
               analysis = x, 
               intercept = analysis$experiment$ineffInterceptActual,
               slope = analysis$experiment$ineffSlopeActual,
               col = actual.col,
               lty = actual.lty
               );
         }
         plotU(analysis = x, xlab = "", ylab = u.ylab);
         if (actual) {
            linesUModel(
               analysis = x,
               umax = analysis$experiment$umax,
               halfsat = analysis$experiment$halfsat,
               col = actual.col,
               lty = actual.lty
               );
         }
         plotSw(analysis = x, xlab = "", ylab = sw.ylab);
         if (actual) {
            linesSwModel(
               analysis = x,
               intercept = analysis$experiment$swInterceptActual,
               slope = analysis$experiment$swSlopeActual,
               col = actual.col,
               lty = actual.lty
               );
         }
         plotTau(analysis = x, xlab = "", ylab = tau.ylab, pch = 1);
         mtext(
            xlab,
            side = 1,
            outer = TRUE,
            line = 1,
            cex = xlab.cex
            );
      }
   );

# HyperbolicAnalysis.plotVf method ####

#' Method plotVf
#' 
#' @exportMethod plotVf

setGeneric(
   name = "plotVf",
   def = function(analysis, ...) { standardGeneric("plotVf") }
   );

#' Plot the uptake velocity analysis
#'
#' @aliases plotVf.HyperbolicAnalysis
#' @exportMethod plotVf

setMethod(
   f = "plotVf",
   signature = "HyperbolicAnalysis",
   definition = function(
      analysis,
      ylim = c(
         min(analysis$metrics$vf),
         max(analysis$metrics$vf, 1 / analysis$vfEstimates$intercept)
         ),
      ylab = bquote(paste(v[f])),
      ...
      ) 
      {
         plot(
            x = analysis$metrics$cefftot, 
            y = analysis$metrics$vf, 
            ylim = ylim,
            ylab = ylab,
            ...
            );
         linesVfModel(analysis);
      }
   );

# HyperbolicAnalysis.linesVfModel method ####

#' Method linesVfModel
#' 
#' @exportMethod linesVfModel

setGeneric(
   name = "linesVfModel",
   def = function(analysis, ...) { standardGeneric("linesVfModel") }
   );

#' Plot lines for the uptake velocity model
#'
#' @aliases linesVfModel.HyperbolicAnalysis
#' @exportMethod linesVfModel

setMethod(
   f = "linesVfModel",
   signature = "HyperbolicAnalysis",
   definition = function(
      analysis,
      intercept = analysis$vfEstimates$intercept,
      slope = analysis$vfEstimates$slope,
      ...
      ) 
      {
         xvals <- seq(
            from = 0,
            to = max(analysis$metrics$cefftot),
            length.out = 30
            )
         ineffModel <- intercept + slope * xvals; 
         lines(
            x = xvals,
            y = 1 / ineffModel,
            ...
            );
      }
   );

# HyperbolicAnalysis.plotIneff method ####

#' Method plotIneff
#' 
#' @exportMethod plotIneff

setGeneric(
   name = "plotIneff",
   def = function(analysis, ...) { standardGeneric("plotIneff") }
   );

#' Plot the inefficiency analysis
#'
#' @aliases plotIneff
#' @exportMethod plotIneff

setMethod(
   f = "plotIneff",
   signature = "HyperbolicAnalysis",
   definition = function(
      analysis,
      xlim = c(
         min(analysis$metrics$cefftot),
         max(analysis$metrics$cefftot)
         ),
      ylim = c(
         min(1 / analysis$metrics$vf),
         max(1 / analysis$metrics$vf, analysis$vfEstimates$intercept)
         ),
      xlab = "Effective concentration",
      ylab = bquote(paste(v[f])),
      pch = 16,
      col = "black",
      actual = TRUE,
      actual.col = "red",
      actual.lty = "dashed",
      ...
      ) 
      {
         ineff <- 1 / analysis$metrics$vf;
         plot(
            x = analysis$metrics$cefftot, 
            y = ineff, 
            xlim = xlim,
            ylim = ylim,
            xlab = xlab,
            ylab = ylab,
            pch = pch,
            col = col,
            ...
            );
         linesIneffModel(analysis);
         if (actual) {
            linesIneffModel(
               analysis = analysis, 
               xlim = xlim,
               intercept = analysis$experiment$vfInterceptActual,
               slope = analysis$experiment$vfSlopeActual,
               col = actual.col,
               lty = actual.lty
               );
         }
      }
   );

# HyperbolicAnalysis.linesIneffModel method ####

#' Method linesIneffModel
#' 
#' @exportMethod linesIneffModel

setGeneric(
   name = "linesIneffModel",
   def = function(analysis, ...) { standardGeneric("linesIneffModel") }
   );

#' Plot lines for the uptake velocity model
#'
#' @aliases linesIneffModel.HyperbolicAnalysis
#' @exportMethod linesIneffModel

setMethod(
   f = "linesIneffModel",
   signature = "HyperbolicAnalysis",
   definition = function(
      analysis,
      xlim = c(
         0,
         max(analysis$metrics$cefftot)
         ),
      grain = 30,
      intercept = analysis$vfEstimates$intercept,
      slope = analysis$vfEstimates$slope,
      col = "black",
      ...
      ) 
      {
         xvals <- seq(
            from = xlim[1],
            to = xlim[2],
            length.out = grain
            )
         ineffModel <- intercept + slope * xvals; 
         lines(
            x = xvals,
            y = ineffModel,
            col = col,
            ...
            );
      }
   );

# HyperbolicAnalysis.plotTau method ####

#' Method plotTau
#' 
#' @exportMethod plotTau

setGeneric(
   name = "plotTau",
   def = function(analysis, ...) { standardGeneric("plotTau") }
   );

#' Plot the tau analysis
#'
#' @aliases plotTau
#' @exportMethod plotTau

setMethod(
   f = "plotTau",
   signature = "HyperbolicAnalysis",
   definition = function(
      analysis,
      xlim = c(
         min(analysis$metrics$cefftot),
         max(analysis$metrics$cefftot)
         ),
      ylim = c(
         min(1 / analysis$metrics$k, analysis$tauEstimates$intercept),
         max(1 / analysis$metrics$k)
         ),
      xlab = "Effective concentration",
      ylab = bquote(paste("\tau")),
      pch = 16,
      col = "black",
      actual = TRUE,
      actual.col = "red",
      actual.lty = "dashed",
      ...
      ) 
      {
         tauNet <- 1 / analysis$metrics$k;
         plot(
            x = analysis$metrics$cefftot, 
            y = tauNet, 
            xlim = xlim,
            ylim = ylim,
            xlab = xlab,
            ylab = ylab,
            pch = pch,
            col = col,
            ...
            );
         linesTauModel(analysis);
         if (actual) {
            linesTauModel(
               analysis = analysis, 
               xlim = xlim,
               intercept = analysis$experiment$tauAmbActual,
               slope = analysis$experiment$tauSlopeActual,
               col = actual.col,
               lty = actual.lty
               );
         }
      }   
   );

# HyperbolicAnalysis.linesTauModel method ####

#' Method linesTauModel
#' 
#' @exportMethod linesTauModel

setGeneric(
   name = "linesTauModel",
   def = function(analysis, ...) { standardGeneric("linesTauModel") }
   );

#' Plot lines for the uptake velocity model
#'
#' @aliases linesTauModel.HyperbolicAnalysis
#' @exportMethod linesTauModel

setMethod(
   f = "linesTauModel",
   signature = "HyperbolicAnalysis",
   definition = function(
      analysis,
      xlim = c(
         0,
         max(analysis$metrics$cefftot)
         ),
      intercept = analysis$tauEstimates$intercept,
      slope = analysis$tauEstimates$slope,
      col = "black",
      ...
      ) 
      {
         tauModel <- intercept + slope * xlim; 
         lines(
            x = xlim,
            y = tauModel,
            col = col,
            ...
            );
      }
   );

# HyperbolicAnalysis.plotU method ####

#' Method plotU
#' 
#' @exportMethod plotU

setGeneric(
   name = "plotU",
   def = function(analysis, ...) { standardGeneric("plotU") }
   );

#' Plot the uptake nonlinear analysis
#'
#' @aliases plotU.HyperbolicAnalysis
#' @exportMethod plotU

setMethod(
   f = "plotU",
   signature = "HyperbolicAnalysis",
   definition = function(
      analysis, 
      ylim = c(
         0,
         max(analysis$uEstimates$umax)
         ),
      ylab = "U",
      ...
      ) 
      {
         plot(
            x = analysis$metrics$cefftot, 
            y = analysis$uEstimates$uamb + analysis$metrics$u, 
            ylim = ylim,
            ylab = ylab,
            ...
            );
         linesUModel(analysis);
      }
   );

# HyperbolicAnalysis.linesUModel method ####

#' Method linesUModel
#' 
#' @exportMethod linesUModel

setGeneric(
   name = "linesUModel",
   def = function(analysis, ...) { standardGeneric("linesUModel") }
   );

#' Plot the line for an uptake kinetic model
#'
#' @aliases linesUModel.HyperbolicAnalysis
#' @exportMethod linesUModel

setMethod(
   f = "linesUModel",
   signature = "HyperbolicAnalysis",
   definition = function(
      analysis,
      umax = analysis$uEstimates$umax,
      halfsat = analysis$uEstimates$halfsat,
      ...
      ) 
      {
         xvals <- seq(
            from = 0,
            to = max(analysis$metrics$cefftot),
            length.out = 30
            )
         lines(
            x = xvals,
            y = hyperbolic(
               umax = umax, 
               halfsat = halfsat, 
               conc = xvals
               ),
            ...
            );
      }
   );

# HyperbolicAnalysis.plotSw method ####

#' Method plotSw
#' 
#' @exportMethod plotSw

setGeneric(
   name = "plotSw",
   def = function(analysis, ...) { standardGeneric("plotSw") }
   );

#' Plot the uptake length analysis
#'
#' @aliases plotSw.HyperbolicAnalysis
#' @exportMethod plotSw

setMethod(
   f = "plotSw",
   signature = "HyperbolicAnalysis",
   definition = function(
      analysis, 
      ylim = c(
         min(analysis$metrics$sw, analysis$swEstimates$intercept),
         max(analysis$metrics$sw)
         ),
      ylab = bquote(paste(s[w])),
      ...
      ) 
      {
         plot(
            analysis$metrics$cefftot, 
            analysis$metrics$sw,
            ylim = ylim,
            ylab = ylab,
            ...
            );
         linesSwModel(analysis);
      }
   );

# HyperbolicAnalysis.linesSwModel method ####

#' Method linesSwModel
#' 
#' @exportMethod linesSwModel

setGeneric(
   name = "linesSwModel",
   def = function(analysis, ...) { standardGeneric("linesSwModel") }
   );

#' Plot the line for the uptake length kintic model
#'
#' @aliases linesSwModel.HyperbolicAnalysis
#' @exportMethod linesSwModel

setMethod(
   f = "linesSwModel",
   signature = "HyperbolicAnalysis",
   definition = function(
      analysis,
      intercept = analysis$swEstimates$intercept,
      slope = analysis$swEstimates$slope,
      ...
      ) 
      {
         maxceff <- max(analysis$metrics$cefftot);
         lines(
            x = c(0, maxceff),
            y = c(
               intercept,
               intercept + slope * maxceff
               ),
            ...
            );
      }
   );

# Class HyperbolicAnalysisMCR defintion and constructor ####

#' Class HyperbolicAnalysisMCR
#' 
#' Class \code{HyperbolicAnalysisMCR} defines an analysis of whole-stream
#' kinetic behavior based on multilevel constant rate (MCR) 
#' solute uptake experiments
#' 
#' @export HyperbolicAnalysisMCR
#' @exportClass HyperbolicAnalysisMCR

HyperbolicAnalysisMCR <- setClass(
   Class = "HyperbolicAnalysisMCR",
   contains = "HyperbolicAnalysis"
   );

#' Constructor method for HyperbolicAnalysis class.
#' 
#' @rdname HyperbolicAnalysisMCR-class
#' @return The object of class \code{HyperbolicAnalysisMCR} created
#'      by the constructor

setMethod(
   f = "initialize",
   signature = "HyperbolicAnalysisMCR",
   definition = function(
      .Object, 
      simulation, 
      analysisWindow,
      conserveColumn = length(simulation$conserveSolute),
      activeColumn = length(simulation$activeSolute),
      reachLengths = NA,
      useRegression = TRUE
      ) 
      {
         numColumns <- length(conserveColumn);
         isSingleColumn <- numColumns == 1;
         indices <- analysisWindow / simulation$outputTimeStep + 1;
         
         if (!(isSingleColumn || useRegression)) {
            time <- rep(analysisWindow, times = numColumns);
         } else {
            time <- analysisWindow;
         }
         
         analysis <- callNextMethod(
            .Object,
            experiment = simulation, 
            time = time,
            indices = indices
            );

         if (isSingleColumn || useRegression)
         {
            if (useRegression)
            {
               consIndex <- conserveColumn[length(conserveColumn)];
               actIndex <- activeColumn[length(activeColumn)];
            }
            conserve <- simulation$conserveSolute[
               indices,
               consIndex
               ];
            active <- simulation$activeSolute[
               indices,
               actIndex
               ];
            
         } else {
            
            conserve <- simulation$conserveSolute[indices, conserveColumn[1]];
            for (i in 2:numColumns) {
               conserve <- c(conserve, simulation$conserveSolute[indices, conserveColumn[i]]);
            }
            active <- simulation$activeSolute[indices, activeColumn[1]];
            for (i in 2:numColumns) {
               active <- c(active, simulation$activeSolute[indices, activeColumn[i]]);
            }
            
         }
         analysis$metrics$conserve <- conserve;
         analysis$metrics$active <- active;
            
         conservebc <- conserve - simulation$conserveBkg;
         analysis$metrics$conservebc <- conservebc;
         
         activebc <- active - simulation$activeBkg;
         analysis$metrics$activebc <- activebc;
         
         activenr <- conservebc * simulation$injectRatio;
         analysis$metrics$activenr <- activenr;
         
         if (isSingleColumn) {
            
            analysis$metrics$travelTime <- simulation$travelTime;
            ceffinject <- sqrt(activebc * activenr);
            
            k <- ( log(simulation$injectRatio)
               - (log(activebc / conservebc)) ) /
               analysis$metrics$travelTime
            
         } else if (useRegression) {
            
            analysis$metrics$travelTime <- simulation$travelTime;
            curvesactive <- as.matrix(
               simulation$activeSolute[indices, activeColumn] - 
                  simulation$activeBkg
               );
            
            ceffinject <- apply(
               curvesactive,
               1,
               function(x) 
                  { 
                  return(exp(sum(log(x)) / length(x)));
                  }
               );
            
            curvescons <- as.matrix(
               simulation$conserveSolute[indices, conserveColumn] - 
                  simulation$conserveBkg
               );
            logcurvesratio <- log(curvesactive / curvescons);
            traveltimes <- reachLengths / simulation$streamVel;
            k <- apply(
               logcurvesratio,
               1,
               function(x) 
                  {
                  lmr <- lm(x ~ traveltimes);
                  return(-lmr$coefficients["traveltimes"])
                  }
               );
            
         } else {
            
            prevIndices <- 1:numColumns;
            
            travelTime <- rep(
               reachLengths[1] / simulation$streamVel, 
               times = numColumns
               );
            ceffinject <- sqrt(activebc[prevIndices] * activenr[prevIndices]);
            k <- ( log(simulation$injectRatio)
               - (log(activebc[prevIndices] / conservebc[prevIndices])) ) /
               travelTime
            
            for (i in 1:(numColumns - 1)) {
               nextIndices <- (i * numColumns + 1):((i + 1) * numColumns);
               
               travelTime <- c(
                  travelTime,
                  rep(
                     reachLengths[i + 1] / simulation$streamVel, 
                     times = numColumns
                     )
                  );
               ceffinject <- c(
                  ceffinject,
                  sqrt(activebc[nextIndices] * activebc[prevIndices])
                  );
               
               k <- c(
                  k,
                  ( (log(activebc[prevIndices] / conservebc[prevIndices]))
                     - (log(activebc[nextIndices] / conservebc[nextIndices])) ) /
                     (travelTime[nextIndices])
                  );
               prevIndices <- nextIndices;
            }
            analysis$metrics$travelTime <- travelTime;
            
         }
         analysis$metrics$ceffinject <- ceffinject;
         analysis$metrics$k <- k;
         
         calcMetrics(analysis);
         
         return(analysis);
      }
   );

# Class HyperbolicAnalysisSlug definition and constructor ####

#' Class HyperbolicAnalysisSlug
#' 
#' Class \code{HyperbolicAnalysisSlug} defines an analysis of whole-stream
#' kinetic behavior based on instantaneous release (slug) 
#' solute uptake experiments
#' 
#' @export HyperbolicAnalysisSlug
#' @exportClass HyperbolicAnalysisSlug

HyperbolicAnalysisSlug <- setClass(
   Class = "HyperbolicAnalysisSlug",
   contains = "HyperbolicAnalysis",
   );

#' Constructor method for HyperbolicAnalysisSlug class.
#' 
#' @rdname HyperbolicAnalysisSlug-class
#' @return The object of class \code{HyperbolicAnalysisSlug} created
#'      by the constructor

setMethod(
   f = "initialize",
   signature = "HyperbolicAnalysisSlug",
   definition = function(
      .Object, 
      simulation, 
      analysisWindow,
      conserveColumn = length(simulation$conserveSolute),
      activeColumn = length(simulation$activeSolute)
      ) 
      {
         startIndex <- trunc(analysisWindow[1] / simulation$outputTimeStep) + 1;
         endIndex <- trunc(analysisWindow[2] / simulation$outputTimeStep);
         indices <- startIndex:endIndex;
         
         time <- simulation$conserveSolute$Time[indices];
         
         analysis <- callNextMethod(
            .Object,
            experiment = simulation, 
            time = time,
            indices = indices
            );

         analysis$releaseTime <- simulation$releaseTime;

         conserve <- simulation$conserveSolute[indices, conserveColumn];
         analysis$metrics$conserve <- conserve;
         
         active <- simulation$activeSolute[indices, activeColumn];
         analysis$metrics$active <- active;
      
         conservebc <- conserve - simulation$conserveBkg;
         analysis$metrics$conservebc <- conservebc;
         
         activebc <- active - simulation$activeBkg;
         analysis$metrics$activebc <- activebc;
         
         activenr <- conservebc * simulation$injectRatio;
         analysis$metrics$activenr <- activenr;
         
         return(analysis);
      }
   );

# Class HyperbolicAnalysisTASCC defintion and constructor ####

#' Class HyperbolicAnalysisTASCC
#' 
#' Class \code{HyperbolicAnalysisTASCC} defines an analysis of whole-stream
#' kinetic behavior based on instantaneous release (slug) 
#' solute uptake experiments and a TASCC analysis 
#' (Covino et al. 2010, Limn. and Ocean.: Methods)
#' 
#' @export HyperbolicAnalysisTASCC
#' @exportClass HyperbolicAnalysisTASCC

HyperbolicAnalysisTASCC <- setClass(
   Class = "HyperbolicAnalysisTASCC",
   contains = "HyperbolicAnalysisSlug",
   );

#' Constructor method for HyperbolicAnalysisTASCC class.
#' 
#' @rdname HyperbolicAnalysisTASCC-class
#' @return The object of class \code{HyperbolicAnalysisTASCC} created
#'      by the constructor

setMethod(
   f = "initialize",
   signature = "HyperbolicAnalysisTASCC",
   definition = function(
      .Object, 
      simulation, 
      analysisWindow,
      ...
      ) 
      {
         analysis <- callNextMethod(
            .Object,
            simulation = simulation, 
            analysisWindow = analysisWindow,
            ...
            );
         
         analysis$metrics$travelTime <- 
            analysis$metrics$time - analysis$releaseTime;
         analysis$metrics$ceffinject <-
            sqrt(analysis$metrics$activebc * analysis$metrics$activenr);
         analysis$metrics$k <- 
            ( log(simulation$injectRatio) 
               - log(analysis$metrics$activebc / analysis$metrics$conservebc) ) /
               (analysis$metrics$travelTime);
         
         length <- length(analysis$metrics$time);
         analysis$metricsOriginal <- data.frame(
            swo = numeric(length = length),
            vfo = numeric(length = length),
            uo = numeric(length = length),
            stringsAsFactors = FALSE
            );
         
         calcMetrics(analysis);
         
         return(analysis);
      }
   );

# Class HyperbolicAnalysisTASCCLagrange defintion and constructor ####

#' Class HyperbolicAnalysisTASCCLagrange
#' 
#' Class \code{HyperbolicAnalysisTASCCLagrange} defines an analysis of whole-stream
#' kinetic behavior based on instantaneous release (slug) 
#' solute uptake experiments and a lagrangian adaptation of a TASCC analysis 
#' 
#' @export HyperbolicAnalysisTASCCLagrange
#' @exportClass HyperbolicAnalysisTASCCLagrange

HyperbolicAnalysisTASCCLagrange <- setClass(
   Class = "HyperbolicAnalysisTASCCLagrange",
   contains = "HyperbolicAnalysisSlug",
   slots = c(
      indices = "numeric"
      )
   );

#' Constructor method for HyperbolicAnalysisTASCCLagrange class.
#' 
#' @rdname HyperbolicAnalysisTASCCLagrange-class
#' @return The object of class \code{HyperbolicAnalysisTASCCLagrange} created
#'      by the constructor

setMethod(
   f = "initialize",
   signature = "HyperbolicAnalysisTASCCLagrange",
   definition = function(
      .Object, 
      simulation, 
      analysisWindow,
      ...
      ) 
      {
         analysis <- callNextMethod(
            .Object,
            simulation = simulation, 
            analysisWindow = analysisWindow,
            ...
            );
         analysis$upstream <- data.frame(
            time = numeric(length(analysis$metrics$time)),
            activebc = numeric(length(analysis$metrics$time)),
            conservebc = numeric(length(analysis$metrics$time))
            );
            
         for (i in 1:length(analysis$metrics$time))
         {
            activebc <- 
               simulation$paths[[i]]$active - simulation$activeBkg;
            conservebc <- 
               simulation$paths[[i]]$conserve - simulation$conserveBkg;
            analysis$upstream$time[i] <- simulation$paths[[i]]$time[1];
            analysis$upstream$activebc[i] <- activebc[1];
            analysis$upstream$conservebc[i] <- conservebc[1];
         }
         analysis$metrics$travelTime <- 
            analysis$metrics$time - analysis$upstream$time;
         analysis$metrics$k <- 
            (log(analysis$upstream$activebc / analysis$upstream$conservebc)
            - log(analysis$metrics$activebc / analysis$metrics$conservebc)) /
               (analysis$metrics$travelTime);
         analysis$metrics$ceffinject <- 
            (analysis$upstream$activebc * analysis$metrics$activebc)^0.5;
         
         calcMetrics(analysis);
         
         return(analysis);
      }
   );
