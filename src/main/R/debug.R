rm(list=ls());
detach("package:strtracer");
remove.packages("strtracer");

install.packages(
   "C:/mount/boxsync/repositories/git/strtracer/src/main/R/strtracer", 
   repos=NULL, 
   type="source"
   );
library(strtracer);
