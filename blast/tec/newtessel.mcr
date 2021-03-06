#!MC 1410
$!VarSet |MFBD| = './'
$!DRAWGRAPHICS FALSE
$!READDATASET  '"|MFBD|/all.dat" '
  READDATAOPTION = NEW
  RESETSTYLE = YES
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
# VARNAMELIST = '"V1" "V2" "V3"'
  VARNAMELIST = '"V1" "V2" "V3" "V4"'
$!PLOTTYPE = CARTESIAN2D
$!TRIANGULATE 
  SOURCEZONES =  [1]
  BOUNDARYZONES =  []
  USEBOUNDARY = NO
  INCLUDEBOUNDARYPTS = NO
  TRIANGLEKEEPFACTOR = 0.25
$!WRITEDATASET  "|MFBD|/triangulatedData.dat"
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  ZONELIST =  [2]
  BINARY = NO
  USEPOINTFORMAT = YES
  PRECISION = 9
  TECPLOTVERSIONTOWRITE = TECPLOTCURRENT
$!RemoveVar |MFBD|
