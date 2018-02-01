setOldClass(Classes = c("Simulation", "Experiment"));

# Class HyperbolicAnalysis definition and constructor ####
HyperbolicAnalysis <- setClass(
   Class = "HyperbolicAnalysis",
   slots = c(
      experiment = "Experiment",
      metrics = "data.frame",
      results = "environment"
      )
   );

setMethod(
   f = "initialize",
   signature = "HyperbolicAnalysis",
   definition = function(
      .Object, 
      experiment, 
      time,
      conserve,
      active,
      conservebc,
      activebc,
      activenr,
      ceffinject,
      k
      ) 
      {
         metricsLength = length(time);
         metrics = new(
            Class = "data.frame",
            data.frame(
               time = time,
               conserve = conserve,
               active = active,
               conservebc = conservebc,
               activebc = activebc,
               activenr = activenr,
               ceffinject = ceffinject,
               cefftot = numeric(length = metricsLength),
               k = k,
               sw = numeric(length = metricsLength),
               vf = numeric(length = metricsLength),
               u = numeric(length = metricsLength),
               stringsAsFactors = FALSE
               )
            );
         metrics$cefftot <- experiment$activeBkg + metrics$ceffinject;
         metrics$sw <- experiment$streamVel / metrics$k;
         metrics$vf <- metrics$k * experiment$streamDepth;
         metrics$u <- metrics$vf * metrics$ceffinject;

         callNextMethod(
            .Object, 
            experiment = experiment, 
            metrics = metrics,
            results = new.env()
            );
      }
   );

# HyperbolicAnalysis.runAll method ####

setGeneric(
   name = "runAll",
   def = function(analysis, ...) { standardGeneric("runAll") }
   );

setMethod(
   f = "runAll",
   signature = "HyperbolicAnalysis",
   definition = function(analysis, ...) {
      
      runSw(analysis, ...);
      runU(analysis, ...);
      runVf(analysis, ...);
   
      }
   
   );

# HyperbolicAnalysis.runSw method ####
setGeneric(
   name = "runSw",
   def = function(analysis) { standardGeneric("runSw") }
   );

setMethod(
   f = "runSw",
   signature = "HyperbolicAnalysis",
   definition = function(analysis) {
      
      # Regression of added solute uptake length vs. concentration
      lmresults <- lm(
        sw ~ cefftot,
        data = analysis@metrics
      );
      env <- analysis@results;
      env$swEstimates <- list(
         intercept = as.numeric(lmresults$coefficients["(Intercept)"]),
         slope = as.numeric(lmresults$coefficients["cefftot"]),
         swamb = numeric(length = 1),
         umax = numeric(length = 1),
         halfsat = numeric(length = 1),
         uamb = numeric(length = 1)
         );
      env$swEstimates$swamb <- env$swEstimates$intercept;
      env$swEstimates$halfsat <- 
         ( (env$swEstimates$intercept + env$swEstimates$slope * analysis@experiment$activeBkg) 
         / env$swEstimates$slope ) - analysis@experiment$activeBkg;
      env$swEstimates$umax <- ( analysis@experiment$discharge * 
         (env$swEstimates$halfsat + analysis@experiment$activeBkg) ) /
         (analysis@experiment$streamWidth * env$swEstimates$slope * env$swEstimates$halfsat)
      env$swEstimates$uamb <- (analysis@experiment$discharge * analysis@experiment$activeBkg) /
            (analysis@experiment$streamWidth * env$swEstimates$swamb);

      }
   );

# HyperbolicAnalysis.runU method ####
setGeneric(
   name = "runU",
   def = function(analysis, ...) { standardGeneric("runU") }
   );

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
               umaxp = analysis@results$swEstimates$umax,
               halfsatp = analysis@results$swEstimates$halfsat
               );
         }
         else
         {
            start = initialEstimates;
         }
         activeBkg <- analysis@experiment$activeBkg;
         
         if (!is.null(fixedParameters$umaxp)) umaxp <- fixedParameters$umaxp;
         if (!is.null(fixedParameters$halfsatp)) halfsatp <- fixedParameters$halfsatp;
         nlsresults <- nls(
            u ~ hyperbolicnet(
              umax = umaxp, 
              halfsat = halfsatp, 
              concadd = ceffinject, 
              concbkg = activeBkg
              ),
            data = analysis@metrics,
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
         env <- analysis@results;
         env$uEstimates <- list(
            umax = umaxest,
            halfsat = halfsatest,
            uamb = hyperbolic(
               umax = umaxest, 
               halfsat = halfsatest, 
               conc = analysis@experiment$activeBkg
               )
            );
      }
   );

# HyperbolicAnalysis.runVf method ####
setGeneric(
   name = "runVf",
   def = function(analysis) { standardGeneric("runVf") }
   );

setMethod(
   f = "runVf",
   signature = "HyperbolicAnalysis",
   definition = function(analysis) {
      ineff <- 1 / analysis@metrics$vf;
      lmresults <- lm(
         ineff ~ cefftot,
         data = analysis@metrics
      );
      intercept = as.numeric(lmresults$coefficients["(Intercept)"]);
      slope = as.numeric(lmresults$coefficients["cefftot"]);
      halfsat = intercept / slope;
      vfambest = 1 / intercept;
      analysis@results$vfEstimates = list(
         intercept = intercept,
         slope = slope,
         vfamb = vfambest,
         umax = (halfsat + analysis@experiment$activeBkg) /
            intercept,
         halfsat = halfsat,
         uamb = vfambest * analysis@experiment$activeBkg
      );
      }
   );

# HyperbolicAnalysis.plot method ####

setGeneric(name = "plot");

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
               intercept = analysis@experiment$vfInterceptActual,
               slope = analysis@experiment$vfSlopeActual,
               col = actual.col,
               lty = actual.lty
               );
         }
         plotU(analysis = x, xlab = "", ylab = u.ylab);
         if (actual) {
            linesUModel(
               analysis = x,
               umax = analysis@experiment$umax,
               halfsat = analysis@experiment$halfsat,
               col = actual.col,
               lty = actual.lty
               );
         }
         plotSw(analysis = x, xlab = "", ylab = sw.ylab);
         if (actual) {
            linesSwModel(
               analysis = x,
               intercept = analysis@experiment$swInterceptActual,
               slope = analysis@experiment$swSlopeActual,
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

setGeneric(
   name = "plotVf",
   def = function(analysis, ...) { standardGeneric("plotVf") }
   );

setMethod(
   f = "plotVf",
   signature = "HyperbolicAnalysis",
   definition = function(
      analysis,
      ylim = c(
         min(analysis@metrics$vf),
         max(analysis@metrics$vf, 1 / analysis@results$vfEstimates$intercept)
         ),
      ylab = bquote(paste(v[f])),
      ...
      ) 
      {
         plot(
            x = analysis@metrics$cefftot, 
            y = analysis@metrics$vf, 
            ylim = ylim,
            ylab = ylab,
            ...
            );
         linesVfModel(analysis);
      }
   );

# HyperbolicAnalysis.linesVfModel method ####

setGeneric(
   name = "linesVfModel",
   def = function(analysis, ...) { standardGeneric("linesVfModel") }
   );

setMethod(
   f = "linesVfModel",
   signature = "HyperbolicAnalysis",
   definition = function(
      analysis,
      intercept = analysis@results$vfEstimates$intercept,
      slope = analysis@results$vfEstimates$slope,
      ...
      ) 
      {
         xvals <- seq(
            from = 0,
            to = max(analysis@metrics$cefftot),
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

setGeneric(
   name = "plotU",
   def = function(analysis, ...) { standardGeneric("plotU") }
   );

setMethod(
   f = "plotU",
   signature = "HyperbolicAnalysis",
   definition = function(
      analysis, 
      ylim = c(
         0,
         max(analysis@results$uEstimates$umax)
         ),
      ylab = "U",
      ...
      ) 
      {
         plot(
            x = analysis@metrics$cefftot, 
            y = analysis@results$uEstimates$uamb + analysis@metrics$u, 
            ylim = ylim,
            ylab = ylab,
            ...
            );
         linesUModel(analysis);
      }
   );

# HyperbolicAnalysis.linesUModel method ####

setGeneric(
   name = "linesUModel",
   def = function(analysis, ...) { standardGeneric("linesUModel") }
   );

setMethod(
   f = "linesUModel",
   signature = "HyperbolicAnalysis",
   definition = function(
      analysis,
      umax = analysis@results$uEstimates$umax,
      halfsat = analysis@results$uEstimates$halfsat,
      ...
      ) 
      {
         xvals <- seq(
            from = 0,
            to = max(analysis@metrics$cefftot),
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

setGeneric(
   name = "plotSw",
   def = function(analysis, ...) { standardGeneric("plotSw") }
   );

setMethod(
   f = "plotSw",
   signature = "HyperbolicAnalysis",
   definition = function(
      analysis, 
      ylim = c(
         min(analysis@metrics$sw, analysis@results$swEstimates$intercept),
         max(analysis@metrics$sw)
         ),
      ylab = bquote(paste(s[w])),
      ...
      ) 
      {
         plot(
            analysis@metrics$cefftot, 
            analysis@metrics$sw,
            ylim = ylim,
            ylab = ylab,
            ...
            );
         linesSwModel(analysis);
      }
   );

# HyperbolicAnalysis.linesSwModel method ####

setGeneric(
   name = "linesSwModel",
   def = function(analysis, ...) { standardGeneric("linesSwModel") }
   );

setMethod(
   f = "linesSwModel",
   signature = "HyperbolicAnalysis",
   definition = function(
      analysis,
      intercept = analysis@results$swEstimates$intercept,
      slope = analysis@results$swEstimates$slope,
      ...
      ) 
      {
         maxceff <- max(analysis@metrics$cefftot);
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

HyperbolicAnalysisMCR <- setClass(
   Class = "HyperbolicAnalysisMCR",
   contains = "HyperbolicAnalysis",
   slots = c(
      indices = "numeric"
      )
   );

setMethod(
   f = "initialize",
   signature = "HyperbolicAnalysisMCR",
   definition = function(
      .Object, 
      simulation, 
      analysisWindow,
      conserveColumn = length(simulation$conserveSolute),
      activeColumn = length(simulation$activeSolute),
      reachLengths = NA
      ) 
      {
         numColumns <- length(conserveColumn);
         isSingleColumn <- numColumns == 1;
         .Object@indices <- analysisWindow / simulation$outputTimeStep + 1;
         
         if (!isSingleColumn) {
            time <- numeric(length = length(analysisWindow) * numColumns)
         }
         time <- analysisWindow;

         if (isSingleColumn)
         {
            conserve <- simulation$conserveSolute[
               .Object@indices,
               conserveColumn
               ];
            active <- simulation$activeSolute[
               .Object@indices,
               activeColumn
               ];
         } else {
            conserve <- simulation$conserveSolute[.Object@indices, conserveColumn[1]];
            for (i in 2:numColumns) {
               conserve <- c(conserve, simulation$conserveSolute[.Object@indices, conserveColumn[i]]);
            }
            active <- simulation$activeSolute[.Object@indices, activeColumn[1]];
            for (i in 2:numColumns) {
               active <- c(active, simulation$activeSolute[.Object@indices, activeColumn[i]]);
            }
         }
         
         conservebc = conserve - simulation$conserveBkg;
         activebc = active - simulation$activeBkg;
         activenr = conservebc * simulation$injectRatio
         
         if (isSingleColumn) {
            ceffinject <- sqrt(activebc * activenr);
            
            k <- ( log(simulation$injectRatio)
               - (log(activebc / conservebc)) ) /
               simulation$travelTime
         } else {
            prevIndices <- 1:numColumns;
            
            ceffinject <- sqrt(activebc[prevIndices] * activenr[prevIndices]);
            k <- ( log(simulation$injectRatio)
               - (log(activebc[prevIndices] / conservebc[prevIndices])) ) /
               (reachLengths[1] / simulation$streamVel)
            
            for (i in 1:(numColumns - 1)) {
               indices <- (i * numColumns + 1):((i + 1) * numColumns);
               
               ceffinject <- c(
                  ceffinject,
                  sqrt(activebc[indices] * activebc[prevIndices])
                  );
               
               k <- c(
                  k,
                  ( (log(activebc[prevIndices] / conservebc[prevIndices]))
                     - (log(activebc[indices] / conservebc[indices])) ) /
                     (reachLengths[i + 1] / simulation$streamVel)
                  );
               prevIndices <- indices;
            }
         }
         
         callNextMethod(
            .Object,
            experiment = simulation, 
            time = time,
            conserve = conserve,
            active = active,
            conservebc = conservebc,
            activebc = activebc,
            activenr = activenr,
            ceffinject = ceffinject,
            k = k
            );
      }
   );