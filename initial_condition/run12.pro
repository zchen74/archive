;Section 8, using IDL
Infile='/panfs/roc/groups/8/grifbake/chen3274/create_initial_condition/ct120101.nc';
   ncdf_read,  result, filename=Infile, /all
 
   data = float(result.SPC_CO2)
;use ctm_make_datainfo to generate datainfo%%
Success = CTM_Make_DataInfo( data, $
ThisDataInfo, $
ModelInfo = CTM_TYPE( 'GEOS5_47L'), $
DiagN = 'IJ-AVG-$', $
Tracer = 1 , $
Tau0 = 271728.d0, $
; note plz use tau_cal.pro in the same directory to get the correct Tau0; must match the timing. 
Unit = 'v/v', $
DIM = [72, 46, 47], $
First = [1,1,1] )
CTM_WriteBpch, ThisDataInfo, FileName='/home/grifbake/chen3274/create_initial_condition/CT20120101'
exit

