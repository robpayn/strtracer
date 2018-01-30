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
   where = globalenv(),
   def = function(analysis, ...) { standardGeneric("runAll") }
   );

setMethod(
   f = "runAll",
   signature = "HyperbolicAnalysis",
   definition = function(analysis, ...) {
      
      runSw(analysis, ...);
      runU(analysis, ...);
   
      }
   
   );

# HyperbolicAnalysis.runSw method ####
setGeneric(
   name = "runSw",
   where = globalenv(),
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
   where = globalenv(),
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
         
         if (isSingleColumn) {
            time = analysisWindow;
         } else {
            time = analysisWindow * length(conserveColumn)
         }
         
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
                  sqrt(activebc[indices] * activenr[prevIndices])
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