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
#' @name runU
#' @rdname runU-methods
#' @exportMethod runU

setGeneric(
   name = "runU",
   def = function(analysis, ...) { standardGeneric("runU") }
   );

#' @rdname runU-methods
#' @aliases runU,HyperbolicAnalysis-method

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
#' @name runVf
#' @rdname runVf-methods
#' @exportMethod runVf

setGeneric(
   name = "runVf",
   def = function(analysis) { standardGeneric("runVf") }
   );

#' @rdname runVf-methods
#' @aliases runVf,HyperbolicAnalysis-method

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

# HyperbolicAnalysis.plot method ####

#' Method plot
#' 
#' @name plot
#' @rdname plot-methods
#' @exportMethod plot

setGeneric(name = "plot");

#' @rdname plot-methods
#' @aliases plot,HyperbolicAnalysis-method
#' @details Plots the results from all inference methods

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
      actual = TRUE,
      actual.col = "red",
      actual.lty = "dashed",
      ...
      ) 
      {
         par(
            mar = c(2, 4.5, 1, 1),
            oma = c(2, 0, 0, 0),
            mfrow = c(3, 1),
            ...
            );
         plotVf(analysis = x, xlab = "", ylab = vf.ylab);
         if (actual) {
            linesVfModel(
               analysis = x, 
               intercept = analysis$experiment$vfInterceptActual,
               slope = analysis$experiment$vfSlopeActual,
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
#' @name plotVf
#' @rdname plotVf-methods
#' @exportMethod plotVf

setGeneric(
   name = "plotVf",
   def = function(analysis, ...) { standardGeneric("plotVf") }
   );

#' @rdname plotVf-methods
#' @aliases plotVf,HyperbolicAnalysis-method

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
#' @name linesVfModel
#' @rdname linesVfModel-methods
#' @exportMethod linesVfModel

setGeneric(
   name = "linesVfModel",
   def = function(analysis, ...) { standardGeneric("linesVfModel") }
   );

#' @rdname linesVfModel-methods
#' @aliases linesVfModel,HyperbolicAnalysis-method

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

# HyperbolicAnalysis.plotU method ####

#' Method plotU
#' 
#' @name plotU
#' @rdname plotU-methods
#' @exportMethod plotU

setGeneric(
   name = "plotU",
   def = function(analysis, ...) { standardGeneric("plotU") }
   );

#' @rdname plotU-methods
#' @aliases plotU,HyperbolicAnalysis-method

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
#' @name linesUModel
#' @rdname linesUModel-methods
#' @exportMethod linesUModel

setGeneric(
   name = "linesUModel",
   def = function(analysis, ...) { standardGeneric("linesUModel") }
   );

#' @rdname linesUModel-methods
#' @aliases linesUModel,HyperbolicAnalysis-method

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
#' @name plotSw
#' @rdname plotSw-methods
#' @exportMethod plotSw

setGeneric(
   name = "plotSw",
   def = function(analysis, ...) { standardGeneric("plotSw") }
   );

#' @rdname plotSw-methods
#' @aliases plotSw,HyperbolicAnalysis-method

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
#' @name linesSwModel
#' @rdname linesSwModel-methods
#' @exportMethod linesSwModel

setGeneric(
   name = "linesSwModel",
   def = function(analysis, ...) { standardGeneric("linesSwModel") }
   );

#' @rdname linesSwModel-methods
#' @aliases linesSwModel,HyperbolicAnalysis-method

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
#' @name HyperbolicAnalysisMCR-class
#' @rdname HyperbolicAnalysisMCR-class
#' @export HyperbolicAnalysisMCR
#' @exportClass HyperbolicAnalysisMCR

HyperbolicAnalysisMCR <- setClass(
   Class = "HyperbolicAnalysisMCR",
   contains = "HyperbolicAnalysis"
   );

#' Constructor method for HyperbolicAnalysisMCR class.
#' 
#' @rdname HyperbolicAnalysisMCR-class

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
            time <- numeric(length = length(analysisWindow) * numColumns)
         }
         time <- analysisWindow;
         
         analysis <- callNextMethod(
            .Object,
            experiment = simulation, 
            time = time,
            indices = indices
            );

         if (isSingleColumn || useRegression)
         {
            
            conserve <- simulation$conserveSolute[
               indices,
               conserveColumn
               ];
            active <- simulation$activeSolute[
               indices,
               activeColumn
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
            
            ceffinject <- sqrt(activebc * activenr);
            
            k <- ( log(simulation$injectRatio)
               - (log(activebc / conservebc)) ) /
               simulation$travelTime
            
         } else if (useRegression) {
            
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
            
            ceffinject <- sqrt(activebc[prevIndices] * activenr[prevIndices]);
            k <- ( log(simulation$injectRatio)
               - (log(activebc[prevIndices] / conservebc[prevIndices])) ) /
               (reachLengths[1] / simulation$streamVel)
            
            for (i in 1:(numColumns - 1)) {
               nextIndices <- (i * numColumns + 1):((i + 1) * numColumns);
               
               ceffinject <- c(
                  ceffinject,
                  sqrt(activebc[nextIndices] * activebc[prevIndices])
                  );
               
               k <- c(
                  k,
                  ( (log(activebc[prevIndices] / conservebc[prevIndices]))
                     - (log(activebc[nextIndices] / conservebc[nextIndices])) ) /
                     (reachLengths[i + 1] / simulation$streamVel)
                  );
               prevIndices <- nextIndices;
            }
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
#' @name HyperbolicAnalysisSlug-class
#' @rdname HyperbolicAnalysisSlug-class
#' @export HyperbolicAnalysisSlug
#' @exportClass HyperbolicAnalysisSlug

HyperbolicAnalysisSlug <- setClass(
   Class = "HyperbolicAnalysisSlug",
   contains = "HyperbolicAnalysis",
   );

#' Constructor method for HyperbolicAnalysisSlug class.
#' 
#' @rdname HyperbolicAnalysisSlug-class
#' 

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
#' @name HyperbolicAnalysisTASCC-class
#' @rdname HyperbolicAnalysisTASCC-class
#' @export HyperbolicAnalysisTASCC
#' @exportClass HyperbolicAnalysisTASCC

HyperbolicAnalysisTASCC <- setClass(
   Class = "HyperbolicAnalysisTASCC",
   contains = "HyperbolicAnalysisSlug",
   );

#' Constructor method for HyperbolicAnalysisTASCC class.
#' 
#' @rdname HyperbolicAnalysisTASCC-class
#' 

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
         
         analysis$metrics$ceffinject <-
            sqrt(analysis$metrics$activebc * analysis$metrics$activenr);
         analysis$metrics$k <- 
            ( log(simulation$injectRatio) 
               - log(analysis$metrics$activebc / analysis$metrics$conservebc) ) /
               (analysis$metrics$time - analysis$releaseTime);
         
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
#' @name HyperbolicAnalysisTASCCLagrange-class
#' @rdname HyperbolicAnalysisTASCCLagrange-class
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
         
         for (i in 1:length(analysis$metrics$time))
         {
            activebc <- 
               simulation$paths[[i]]$active - simulation$activeBkg;
            conservebc <- 
               simulation$paths[[i]]$conserve - simulation$conserveBkg;
            logy <- log(activebc / conservebc);
            x <- simulation$paths[[i]]$time;
            lmresults <- lm(logy ~ x);
            analysis$metrics$k[i] <- -lmresults$coefficients["x"];
            analysis$metrics$ceffinject[i] <- exp(mean(log(activebc)));
         }
         
         calcMetrics(analysis);
         
         return(analysis);
      }
   );
