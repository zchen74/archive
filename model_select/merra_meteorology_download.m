cd /home/grifbake/chen3274/shared/data/geos-chem/ExtData/GEOS_4x5
for i=2014:2019;
    ii=num2str(i);
    
    for j=1:12;
        
        if j<10; jj=strcat('0',num2str(j));end;
        if j>=10; jj=num2str(j);end;
        
      
        
        filename=strcat('wget -r -nH --cut-dirs=5 ftp://ftp.as.harvard.edu/gcgrid/geos-chem/data/ExtData/GEOS_4x5/MERRA2/',ii,'/',jj,'/ ./');
        unix(filename);
    end
end
