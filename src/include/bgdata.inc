      integer maxbgmodels
      parameter (maxbgmodels=14) ! Yuanfu Xie added 14 for processing GRAPES
      real missingflag
      parameter (missingflag=1.0e37)
      real esat(15000:45000), es(15000:45000) 	
      common /estab/ esat,es
