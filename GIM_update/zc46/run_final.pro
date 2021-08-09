; Written by Zichong Chen, Dec 24th, 2018
; Johns Hopkins University
 FOR i=1,365 DO BEGIN &$

Infile='/home/grifbake/chen3274/ntei/zc46/zichong'+ strcompress(i) + '.nc' &$
Infile=STRJOIN(STRSPLIT(Infile, /EXTRACT), '') &$

   ncdf_read,  result, filename=Infile, /all &$
 
   data = float(result.CO2SRCE) &$
;use ctm_make_datainfo to generate datainfo%%
Success = CTM_Make_DataInfo( data, $
ThisDataInfo, $
ModelInfo = CTM_TYPE( 'GEOS5_47L'), $
DiagN = 'CO2-SRCE', $
Tracer = 3 , $
Tau0 = 131448.d0+i*24.d0, $
Unit = 'molec/cm2/s', $
DIM = [72, 46, 1], $
First = [1,1,1] ) &$


ii=STRJOIN(STRSPLIT(strcompress(i), /EXTRACT), '') &$
IF (i LT 10) THEN BEGIN &$
ii='00' + ii &$
ENDIF ELSE IF (i LT 100) THEN BEGIN &$
ii='0' + ii &$
ENDIF &$
CTM_WriteBpch, ThisDataInfo, FileName='/home/grifbake/chen3274/ntei/zc46/CO2.daily.geos.4x5.' + ii &$
ENDFOR
exit
