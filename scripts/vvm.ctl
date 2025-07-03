DSET ^0.nc
DTYPE netcdf
OPTIONS template
TITLE NetCDF Data for GrADS
UNDEF -9999.0
XDEF 502 LINEAR 0 1   # Update if x-spacing isn’t uniform
YDEF 1 LINEAR 1 1     # Single y-point for 2D model
ZDEF 32 LEVELS 0., 44.7368, 165.789, 328.947, 534.211, 781.579, 1071.05, 1402.63, 1776.32, 2192.11, 2650., 3150., 3692.11, 4276.32, 4902.63, 5571.05, 6281.58, 7034.21, 7828.95, 8665.79, 9544.74, 10465.8, 11428.9, 12434.2, 13481.6, 14571.1, 15702.6, 16876.3, 18092.1, 19350., 20650., 21000.
TDEF 6000 LINEAR 00Z31dec1999 2mn # date

VARS 13
th=>th 32 t,x,z theta 
u=>u 32 t,x,z u 
w=>w 32 t,x,z w 
zeta=>zeta 32 t,x,z zeta 
qv=>qv 32 t,x,z qv 
qc=>qc 32 t,x,z qc
qr=>qr 32 t,x,z qr 
qitot=>qitot 32 t,x,z qitot 
precip=>precip 1 t,x precip
ubarTop=>ubarTop 1 t ubarTop
radiation_heating_rate=>rhr 32 t,x,z rhr
heatflux_sfc=>htflx 1 t,x precip
waterflux_sfc=>wtflx 1 t,x precip
ENDVARS

