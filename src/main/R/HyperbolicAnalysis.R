HyperbolicAnalysis <- function(
   experiment, 
   metricsLength
   )
{
   analysis <- new.env();

   analysis$metrics <- data.frame(
      time = numeric(length = metricsLength),
      conserve = numeric(length = metricsLength),
      active = numeric(length = metricsLength),
      conservebc = numeric(length = metricsLength),
      activebc = numeric(length = metricsLength),
      ceffinject = numeric(length = metricsLength),
      cefftot = numeric(length = metricsLength),
      activenr = numeric(length = metricsLength),
      k = numeric(length = metricsLength),
      sw = numeric(length = metricsLength),
      vf = numeric(length = metricsLength),
      u = numeric(length = metricsLength),
      stringsAsFactors = FALSE
      );

   analysis$experiment <- experiment;
   
   class(analysis) <- c("HyperbolicAnalysis", class(analysis));
   return(analysis);   
}

HyperbolicAnalysisTASCCField <- function(experiment)
{
   length <- length(experiment$conserveSolute$Time);
   analysis <- HyperbolicAnalysis(
      experiment,
      length
      );
   
   analysis$metrics$time <- experiment$conserveSolute$Time;
   analysis$metrics$conserve <- experiment$conserveSolute$conserve;
   analysis$metrics$active <- experiment$activeSolute$active;
   
   analysis$metrics$conservebc <- 
      analysis$metrics$conserve - experiment$conserveBkg;
   analysis$metrics$activebc <- 
      analysis$metrics$active - experiment$activeBkg;
   analysis$metrics$activenr <- 
      analysis$metrics$conservebc * experiment$injectRatio;
   
   analysis$metricsOriginal <- data.frame(
      swo = numeric(length = length),
      vfo = numeric(length = length),
      uo = numeric(length = length),
      stringsAsFactors = FALSE
      );
   
   analysis$metrics$ceffinject <-
      sqrt(analysis$metrics$activebc * analysis$metrics$activenr);
   analysis$metrics$k <- 
      ( log(experiment$injectRatio) 
         - log(analysis$metrics$activebc / analysis$metrics$conservebc) ) /
         (analysis$metrics$time - experiment$releaseTime);
   
   calcMetrics(analysis);
   
   return(analysis);
}

HyperbolicAnalysisMultilevel <- function(
   simulation,
   analysisWindow,
   conserveColumn = length(simulation$conserveSolute),
   activeColumn = length(simulation$activeSolute)
   )
{
   analysis <- HyperbolicAnalysis(
      experiment = simulation, 
      metricsLength = length(analysisWindow)
      );
   
   analysis$metrics$time <- analysisWindow;
   analysis$indeces <- 
      analysis$metrics$time / simulation$outputTimeStep + 1;
   analysis$metrics$conserve <- simulation$conserveSolute[
      analysis$indeces,
      conserveColumn
      ];
   analysis$metrics$active <- simulation$activeSolute[
      analysis$indeces,
      activeColumn
      ];
   
   analysis$metrics$conservebc <- 
      analysis$metrics$conserve - simulation$conserveBkg;
   analysis$metrics$activebc <- 
      analysis$metrics$active - simulation$activeBkg;
   analysis$metrics$activenr <- 
      analysis$metrics$conservebc * simulation$injectRatio;
   
   analysis$metrics$ceffinject <-
      sqrt(analysis$metrics$activebc * analysis$metrics$activenr);
   analysis$metrics$k <-
      ( log(simulation$injectRatio)
         - (log(analysis$metrics$activebc / analysis$metrics$conservebc)) ) /
         simulation$travelTime;
   
   calcMetrics(analysis);

   return(analysis);
}

HyperbolicAnalysisSlug <- function(
   simulation,
   analysisWindow,
   conserveColumn = length(simulation$conserveSolute),
   activeColumn = length(simulation$activeSolute)
   )
{
   startIndex <- trunc(analysisWindow[1] / simulation$outputTimeStep) + 1;
   endIndex <- trunc(analysisWindow[2] / simulation$outputTimeStep);

   analysis <- HyperbolicAnalysis(
      experiment = simulation, 
      metricsLength = (endIndex - startIndex) + 1
      );
   
   analysis$startIndex <- startIndex;
   analysis$endIndex <- endIndex;

   analysis$releaseTime <- simulation$releaseTime;

   analysis$metrics$time <- 
      simulation$conserveSolute$Time[analysis$startIndex:analysis$endIndex];

   analysis$metrics$conserve <- simulation$conserveSolute[
      analysis$startIndex:analysis$endIndex,
      conserveColumn
      ];
   analysis$metrics$active <- simulation$activeSolute[
      analysis$startIndex:analysis$endIndex,
      activeColumn
      ];

   analysis$metrics$conservebc <- 
      analysis$metrics$conserve - simulation$conserveBkg;
   analysis$metrics$activebc <- 
      analysis$metrics$active - simulation$activeBkg;
   analysis$metrics$activenr <- 
      analysis$metrics$conservebc * simulation$injectRatio;

   return(analysis);
}

HyperbolicAnalysisTASCC <- function(
   simulation,
   analysisWindow,
   ...
   )
{
   analysis <- HyperbolicAnalysisSlug(simulation, analysisWindow, ...);
   
   length <- length(analysis$metrics[[1]]);
   analysis$metricsOriginal <- data.frame(
      swo = numeric(length = length),
      vfo = numeric(length = length),
      uo = numeric(length = length),
      stringsAsFactors = FALSE
      );
   
   analysis$metrics$ceffinject <-
      sqrt(analysis$metrics$activebc * analysis$metrics$activenr);
   analysis$metrics$k <- 
      ( log(simulation$injectRatio) 
         - log(analysis$metrics$activebc / analysis$metrics$conservebc) ) /
         (analysis$metrics$time - analysis$releaseTime);
   
   calcMetrics(analysis);
  
   return(analysis);
}

HyperbolicAnalysisLagrangeTASCC <- function(
   simulation,
   analysisWindow,
   ...
   )
{
   analysis <- HyperbolicAnalysisSlug(simulation, analysisWindow, ...);
   
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

calcMetrics <- function(analysis, ...)
{
   UseMethod("calcMetrics", analysis);
}

calcMetrics.HyperbolicAnalysis <- function(analysis)
{
   analysis$metrics$cefftot <- analysis$experiment$activeBkg + analysis$metrics$ceffinject;
   analysis$metrics$sw <- analysis$experiment$streamVel / analysis$metrics$k;
   analysis$metrics$vf <- analysis$metrics$k * analysis$experiment$streamDepth;
   analysis$metrics$u <- analysis$metrics$vf * analysis$metrics$ceffinject;
}

run <- function(analysis, ...)
{
   UseMethod("run", analysis);
}

run.HyperbolicAnalysis <- function(
   analysis,
   fixedParameters = list(umaxp = NULL, halfsatp = NULL),
   initialEstimates = NULL
   )
{
   # Regression of added solute uptake length vs. concentration
   lmresults <- lm(
     sw ~ cefftot,
     data = analysis$metrics
   );
   intercept <- as.numeric(lmresults$coefficients["(Intercept)"]);
   slope <- as.numeric(lmresults$coefficients["cefftot"]);
   swambest <- intercept;
   halfsatest <- ((intercept + slope * analysis$experiment$activeBkg) / slope) - 
      analysis$experiment$activeBkg;
   analysis$swEstimates <- list(
      intercept = intercept,
      slope = slope,
      swamb = swambest,
      uamb = (analysis$experiment$discharge * analysis$experiment$activeBkg) /
         (analysis$experiment$streamWidth * swambest),
      halfsat = halfsatest,
      umax = (analysis$experiment$discharge * (halfsatest + analysis$experiment$activeBkg)) /
         (analysis$experiment$streamWidth * slope * halfsatest)
   );

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

   # Linear regression of 1/vf vs. concentration
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
      umax = (halfsat + analysis$experiment$activeBkg) /
         intercept,
      halfsat = halfsat,
      vfamb = vfambest,
      uamb = vfambest * analysis$experiment$activeBkg
   );
}

plot.HyperbolicAnalysis <- function(
   analysis,
   device = "default",
   width = 8,
   height = 6,   
   xfactor = 1,
   yfactor = 1,
   xlim = c(
      0,
      max(analysis$metrics$cefftot)
      ),
   ylim = c(
      0,
      max(analysis$metrics$u)
      ),
   xlab = "Concentration",
   ylab = "Net uptake",
   ...
   )
{
   createDevice(device, width, height);
   par(...);
   createBlankPlot(
      xlim = xlim * xfactor, 
      ylim = ylim * yfactor, 
      xlab = xlab, 
      ylab = ylab
      );
   points(
      x = analysis$metrics$cefftot * xfactor, 
      y = analysis$metrics$u * yfactor
      );
}

plotUptakeEstimate <- function(analysis, ...)
{
   UseMethod("plotUptakeEstimate", analysis);
}

plotUptakeEstimate.HyperbolicAnalysis <- function(
   analysis, 
   device = "default",
   width = 8,
   height = 6,   
   xfactor = 1,
   yfactor = 1,
   xlim = c(
      0,
      max(analysis$metrics$cefftot)
      ),
   ylim = c(
      0,
      max(
         analysis$metrics$u + analysis$uEstimates$uamb,
         if (length(actualModel) == 2) 
            hyperbolic(
               actualModel["umax"], 
               actualModel["halfsat"], 
               max(analysis$metrics$cefftot)
               ) 
         else 0
         ) 
      ),
   xlab = "Concentration",
   ylab = "Uptake",
   actualModel = numeric(length = 0),
   ...
   )
{
   createDevice(device, width, height);
   par(...);
   createBlankPlot(
      xlim = xlim * xfactor, 
      ylim = ylim * yfactor, 
      xlab = xlab, 
      ylab = ylab
      );
   points(
      x = analysis$metrics$cefftot * xfactor, 
      y = (analysis$uEstimates$uamb + analysis$metrics$u) * yfactor,
      pch = 16,
      col = "red"
      );
   xvals <- seq(
      from = 0, 
      max(analysis$metrics$cefftot),
      length.out = 30
      );
   lines(
      x = xvals * xfactor,
      y = hyperbolic(
         umax = analysis$uEstimates$umax, 
         halfsat = analysis$uEstimates$halfsat, 
         conc = xvals
         ) * yfactor
      );
   if (length(actualModel) == 2)
   {
      lines(
         x = xvals * xfactor,
         y = hyperbolic(
            umax = actualModel["umax"], 
            halfsat = actualModel["halfsat"], 
            conc = xvals
            ) * yfactor,
         lty = "dashed"
         );
   }
}

plotVfEstimate <- function(analysis, ...)
{
   UseMethod("plotVfEstimate", analysis);
}

plotVfEstimate.HyperbolicAnalysis <- function(
   analysis, 
   device = "default",
   width = 8,
   height = 6,   
   xfactor = 1,
   yfactor = 1,
   xlim = c(
      0,
      max(analysis$metrics$cefftot)
      ),
   ylim = c(
      min(
         1 / max(analysis$metrics$vf),
         analysis$vfEstimates$intercept,
         if (length(actualModel) == 2) 
            actualModel["vfint"]
         ),
      max(
         1 / analysis$metrics$vf,
         if (length(actualModel) == 2) 
            actualModel["vfint"] +
               actualModel["vfslope"] * max(analysis$metrics$cefftot)
         else 0
         ) 
      ),
   xlab = "Concentration",
   ylab = bquote(paste(
      v[f]^-1
      )),
   actualModel = numeric(length = 0),
   ...
   )
{
   createDevice(device, width, height);
   par(...);
   createBlankPlot(
      xlim = xlim * xfactor, 
      ylim = ylim * yfactor, 
      xlab = xlab, 
      ylab = ylab
      );
   points(
      x = analysis$metrics$cefftot * xfactor, 
      y = (1 / analysis$metrics$vf) * yfactor,
      pch = 16,
      col = "red"
      );
   lines(
      x = xlim,
      y = c(
         analysis$vfEstimates$intercept,
         analysis$vfEstimates$intercept +
            analysis$vfEstimates$slope * xlim[2]
         )
      );
   if (length(actualModel) == 2)
   {
      lines(
         x = xlim,
         y = c(
            actualModel["vfint"],
            actualModel["vfint"] +
               actualModel["vfslope"] * xlim[2]
            ),
         lty = "dashed"
         );
   }
}
