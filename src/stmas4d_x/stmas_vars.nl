&STMAS_variables
STMAS_maxobs=200000
STMAS_numvars=6
STMAS_varnames='U3','V3','W3','P3','T3','SH'
STMAS_radius=
100000.0,100000.0,1800.0,100.0,2.0e0,0.1,
100000.0,100000.0,1800.0,100.0,2.0e0,0.1
100000.0,100000.0,1800.0,100.0,2.0e0,0.1
100000.0,100000.0,1800.0,100.0,4.0e4,0.1
100000.0,100000.0,1800.0,100.0,2.0e0,0.1
30000.0,30000.0,1800.0,100.0,2.0,0.1
STMAS_inc=10.0,10.0,10.0,10.0,8.0,5.0
STMAS_thresholds = 10.0,10.0,10.0,10.0,10.0,5.0
STMAS_stddev = 10.0,10.0,10.0,2.5,2.5,4.0
STMAS_penal=0.5,1.,3.e-2,5.e+5,2.e+4,2.e+2
STMAS_smooth=0.0,0.0,0.0,0.0,0.0,0.0

c numvars         Number of analysis variables;
c varnames        The names of these analysis variables, U3, V3, OM, T3, SH, and HT 
c                 are in lga file. PRES can be the verical grid. One can get PRES from the call 
c                 to get_pres_3d, and convert OM (dpres/dt) to w (m/s) by calling the function omega_to_w(omega,pres);
c                 In order to run a multi-variate analysis, users have
c                 to arrange their variable names as UU,VV,WW,PRES,TT,QQ; Hongli Jiang 6/10/2011
c STMASFC_varnames surface pressure. HJ 6/22/2011
c radius          Influence radius of observations in each directions.
c                 Each variable has different radius in x,y,t,z,background
c                 and land-water are considered (meters, seconds);
c inc             increment sizes for controling how obs mapped to grid following bkgd patterns;
c thresholds      Threshold values for rejecting obs that has larger difference
c                 from background than the thresholds.
c stddev          Factors of standard deviations for rejecting observation data stma
c penal           penality coefficient cooresponing to varnames. 
