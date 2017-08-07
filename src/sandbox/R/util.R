comparisonPlot <- function(
   solstore, solotis, time, timeOtis, locations, locationsOtis, colors, types, xlim, ylim, legend
)
{
   par(mar=c(5,5,4,2), xpd=TRUE);
   plot(
      x = 0, 
      xlab = "Time (sec)",
      ylab = expression(paste(
         "Concentration (g  ",
         m^-3,
         ")"
      )),
      xlim = xlim,
      ylim = ylim,
      type = "n"
   );
   
   for (i in 1:length(locations))
   {
      lines(
         x = time, 
         y = solstore[,sprintf("matrix.cell%s.consConc",locations[i])], 
         col = colors[i],
         lty = types["ktrans"]
      );
      lines(
         x = timeOtis,
         y = solotis[,locationsOtis[i]],
         col = colors[i],
         lty = types["otis"]
      );
   }
   
   legend(
      x = xlim[2] * 0.8,
      y = max(
         solstore[,sprintf("matrix.cell%s.consConc",locations[1])],
         solotis[,locationsOtis[1]]
      ),
      bty = "n",
      legend = c(paste(legend,"NEO"),paste(legend,"OTIS")),
      col = colors,
      lty = c(
         rep(types["ktrans"], length(legend)),
         rep(types["otis"], length(legend))
      )
   );
}
