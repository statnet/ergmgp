#  File R/zzz.R in package ergm.userterms, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2012-2019 Statnet Commons
#######################################################################
.onAttach <- function(lib, pkg){
  sm <- statnetStartupMessage("ergmgp", c("statnet"), FALSE)
  if(!is.null(sm)){
    packageStartupMessage(sm)
  }
}



