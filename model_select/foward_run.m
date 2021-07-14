%% run the model, i.e., from X to h(X)%


%load necessary X; it could be env drivers, or could be fluxes or flux component in TRENDY models%


%we convert X to netcdf format first, and then convert to bpch format, the format needed for GC run%
%why do not I directly convert X to bpch directly -- because I dont know how%

load ~/model_select/TRENDY/orchidee18.mat XX;
m=46*72*365;
X=XX;  
ntimes=365;
X(isnan(X))=0;

for ii=1:1,  
for jj=1:7,  %change the biomes each round%, by biome
j=(ii-1)*7+jj;
Qx=X(:,j);

str=num2str(jj);

temp=Qx.*6.02e+13; %molec/cm2/s is the unit needed from GC model%
temp=reshape(temp,46,72,ntimes);
for i=1:ntimes,
    
    ncid = netcdf.open(strcat('/panfs/roc/groups/8/grifbake/chen3274/model_select/mm',str,'/zichong',num2str(i),'.nc'),'WRITE'); %i indicate itasca%
    varid = netcdf.inqVarID(ncid,'CO2SRCE');
    netcdf.putVar(ncid,varid,temp(:,:,i)'); %Read the NOTE in the inversion_run.m script,i.e., transpose here to agree with the format requirement by the emission input netcdf file%
    netcdf.close(ncid);
end
clear temp;


unix(strcat('/soft/idl/8.0/bin/idl /home/grifbake/chen3274/model_select/mm',str,'/run_final.pro')); %convert to bpch file%

%% Step 2. run the GC forward model, to get Hs;
unix(strcat('sbatch /home/grifbake/chen3274/bbyrne/runs/v8-02-01/m',str,'/run.4'))

end
end
