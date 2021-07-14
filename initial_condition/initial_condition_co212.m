%-------------------------------------------------------------------------------------------------------%
% SCRIPT: initial_condition.m 										%
% PURPOSE: Create an initial condition for GEOS-Chem. Use CarbonTracker-CO2 from Sept 1, 2012		%
% as the initial condition. Then put in the format of GEOS-Chem restart files.				%
% Z. Chen, Dec 23rd, 2018										%
% This script is adapted from the Carbon Tracker-Ch4 script by Dr. Scot Miller					
%-------------------------------------------------------------------------------------------------------%
%% Part 1: Get the inf needed from GEOS-Chem model;
filename='/home/grifbake/chen3274/create_initial_condition/ctm.nc'; %Ap Bp are just constant for the 47 layers%
ap=ncread(filename,'Ap');
bp=ncread(filename,'Bp');


gfilename='/home/grifbake/chen3274/shared/data/geos-chem/ExtData/GEOS_4x5/GEOS_FP/2014/09/GEOSFP.20140901.I3.4x5.nc';
sp=ncread(gfilename,'PS',[1 1 1],[72 46 1]); %the unit is hPa!!%
g_lon0=ncread(gfilename,'lon');
g_lat0=ncread(gfilename,'lat');



press_bottom=[];
for i=1:size(ap,1),
press_bottom(:,:,i)=ap(i)+bp(i).*sp; %hPa%
end
clear ap bp sp;
max(max(max(press_bottom))) 
%%must be larger than 1E3 hPa.
%% 
g_press = zeros(72,46,47);
        for L = 1:47;
        g_press(:,:,L) = 0.5 .* press_bottom(:,:,L) + 0.5 .* press_bottom(:,:,L+1);
        end;
        clear press_bottom;	


ctfile='~/create_initial_condition/CT-NRT.v2015-2.molefrac_glb3x2_2012-09-01.nc';

	ct_lon0 = ncread(ctfile,'lon');
	ct_lat0 = ncread(ctfile,'lat');
    co2 = ncread(ctfile,'co2');
	% Pull out the first time period
	co2=co2(:,:,:,1);
	% Read in pressure coordinates
	% Pressure is in units of hPa
	press = ncread(ctfile,'pressure');
	press = press(:,:,:,1)/100;
   
    %% To avoid nan values
    co2(:,91,:)   = co2(:,90,:);
	co2(121,:,:)   = co2(120,:,:);
    co2(:,2:92,:)=co2;
    co2(2:122,:,:)=co2;
    co2(:,:,26)=co2(:,:,25);
    press(:,91,:) = press(:,90,:);
    press(121,:,:) = press(120,:,:);
    
    press(:,2:92,:)=press;
    press(2:122,:,:)=press;
    

        
        %% 
 	ct_lat0       = [-90;ct_lat0; 90];
 	ct_lon0       = [-180;ct_lon0; 180];

	% Convert lats and lons into a grid
	g_lon  = repmat(g_lon0,1,46);
	g_lat  = repmat(g_lat0',72,1);
        ct_lon  = repmat(ct_lon0,1,92); %NOW WE KNOW WHY LOL%
        ct_lat  = repmat(ct_lat0',122,1);

	% Interpolate CT-CH4 to GEOS-Chem levels
	% interp3(X,Y,Z,V,Xq,Yq,Zq)
	% X,Y,Z are the coordinates. They can be in 'meshgrid' format
	% Here's what I'm going to do:
	% 1. Regrid vertical level of CT-CO2 to the GEOS-Chem lat-lon grid
	% 2. Loop over reach lat-lon grid box in CT-CH4 and interpolate to GEOS-Chem vertical levels
	co22   = [];
	press2 = [];
	for j = 1:size(press,3);
	co22(:,:,j) = interp2(ct_lon',ct_lat',co2(:,:,j)',g_lon',g_lat')';
	press2(:,:,j) = interp2(ct_lon',ct_lat',press(:,:,j)',g_lon',g_lat')';
	end;
    
    co22(:,:,2:27)=co22;
    co22(:,:,28)=co22(:,:,27);
     press2(:,:,2:27)=press2;
    press2(:,:,1)=2000;
    press2(:,:,28)=0;
    

	% Now loop over each lat-lon grid point and interpolate the vertical levels
	co23 = [];
	for j = 1:length(g_lon0);
	for k = 1:length(g_lat0);
	co23(j,k,:) = interp1(reshape(press2(j,k,1:27),27,1),reshape(co22(j,k,1:27),27,1),reshape(g_press(j,k,:),47,1));	
	end;
	end;	
	co2 = co23; % Re-name the interpolated output
	clear press2 co22 co23; % Clear unneeded objects

	% Convert to mol/mol
	co2 = co2./1e6;
	co2 = reshape(co2,72,46,47,1);


%------------------------------%
% Write the new restart file   %
%------------------------------%

	disp('Write new restart file');

	% Format of the GEOS-Chem restart files:
	% Restart file naming convention: GEOSChem_restart.201001020000.nc	

	% Components of the GEOS-Chem restart file:
	% 2 variables (excluding dimension variables):
        % double AREA[lon,lat]   
        % long_name: Grid box area
        % units: m2
        % float SPC_CH4[lon,lat,lev,time]   
        % long_name: SPC_CH4
        %  units: mol mol-1
        %   averaging_method: Instantaneous
        %   _FillValue: -9.99999984824321e+30

        % 4 dimensions:
        % lon  Size:144
        %    long_name: Longitude
        %    units: degrees_east
        % lat  Size:91
        %   long_name: Latitude
        %   units: degrees_north
        % lev  Size:47
        %  long_name: GEOS-Chem level
        %   units: unitless
        % time  Size:1
        %  long_name: Time
        %   units: hours since 2010-01-02 00:00:00 GMT

        % 3 global attributes:
        % title: GEOSChem  restart
        % history: NC_CREATE.F90
        % format: netCDF-4	

	inpath='/home/grifbake/chen3274/create_initial_condition/';
    infile  = strcat(inpath,'zcz.nc');
	outfile = strcat(inpath,'ct120901.nc');
	S = ncinfo(infile);
	ncwriteschema(outfile,S);
	S = ncinfo(outfile);
	ncwrite(outfile,'SPC_CO2',co2); %mol/mol%
	ncwrite(outfile,'lon', ncread(infile,'lon'));
	ncwrite(outfile,'lat', ncread(infile,'lat'));
    ncwrite(outfile,'time',ncread(infile,'time'));
    ncwrite(outfile,'lev',ncread(infile,'lev'));
    
    %unix('idl ~/create_initial_condition/run.pro')
    

%-------------------------------------------------------------------------------------------------------%
% END OF SCRIPT												%
%-------------------------------------------------------------------------------------------------------%
