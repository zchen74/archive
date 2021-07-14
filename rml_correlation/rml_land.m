%---------------------------------------------------------------------------------------------------------------%
% SCRIPT: RML_carbontracker.m											%
% PURPOSE: Run kriging on Carbon Tracker to estimate the covariance matrix parameters of the CT estimate.	%
% S. Miller, Aug. 5, 2018											%
%														%
%---------------------------------------------------------------------------------------------------------------%


%-------------------%
% FUNCTION NOTES:   %
%-------------------%



%-----------------------%
% SET REQUIRED INPUTS   %
%-----------------------%

	% Set path to carbon tracker inputs
    load /home/grifbake/chen3274/testout_L13/data/landmap.mat landmap;
     landmap(landmap==2)=0; %0 for ocean and 1 for land, 2 is for ocean%
   
    create_distmat_land   
	distmat = 0.5 .* (distmat + distmat'); % Make sure the matrix is symmetric
	E  = distmat;
	clear distmat;
     
     
	ctpath = '/home/grifbake/chen3274/ffm/CT/3hr/ctdata/';

        % Create empty output object
        thetaall = [];

	% Loop over each month of the year
	ymlist = {'201501','201502','201503','201504','201505','201506','201507','201508','201509','201510','201511','201512'};
        %ymlist = {'201501','201507'};
	for ym = ymlist;

	ym = ym{1};
	thetatemp = [0 0 0];

	% Set the three hour chunk of the day (value can range from 1 to 8)
	for timechunk = 1:1;

	disp('-------------------------');
	disp('Time chunk');
	disp(num2str(timechunk));


%---------------------------------------%
% Process inputs from CarbonTracker     %
%---------------------------------------%

	flux = [];

	% Loop over each day of the month
	ntimes = 30;
	if strcmp(ym,'201502'); ntimes = 28; end;
	for days = 1:ntimes;

	% Read in the CarbonTracker file
	if days<10;
	infile = strcat(ctpath,'CT2017.flux1x1.',num2str(ym),'0',num2str(days),'.nc'); %ym is alreay a string%
	end;

	if days>9;
	infile = strcat(ctpath,'CT2017.flux1x1.',num2str(ym),num2str(days),'.nc');
	end;

	flux1 = ncread(infile,'bio_flux_opt') + ncread(infile,'fossil_flux_imp') + ncread(infile,'fire_flux_imp'); 

	% Only keep fluxes for the designated time of day
	flux1 = nanmean(flux1,3); %daily mean%
    
    flux1=flux1'; %transpose%
    
clon=-179.5:179.5; clon=repmat(clon,180,1);
clat=-89.5:89.5; clat=repmat(clat',1,360);
load E:\DBAKER\GC_LonLat.mat lon lat; glon=lon;glat=lat;
temp=interp2(clon,clat,flux1,glon,glat);  
flux1=temp; clear temp;    
     
	% Only keep the land regions
	flux1 = flux1(landmap==1); %note there will be only just one column%
	
	% Append to fluxes from existing days
	flux = [flux; flux1];  %make sure it's only 1 column%

	end; % End of days loop


%----------------------------------------------------%
% Create the temporal and spatial distance matrices  %
%----------------------------------------------------%

	% Create the time distance matrix; simply create D%
	days = 1:ntimes; %1:30%
	days = days';
	days = days * ones(1,length(days));
	D = abs(days - days');
	clear days;


%-----------------------------------------------------------------%
% Convert the unit on flux to make the problem more well-behaved  %
%-----------------------------------------------------------------%

	% Convert from grams to micrograms %be careful here zc%
	flux = 1e6 .* flux; %umol/m2/s%


%---------------%
% Create X      %
%---------------%

	X = ones(length(flux),1);


%---------------------------%
% Launch the RML script     %
%---------------------------%

	theta0 = [std(flux),8,1800];

        theta0 = [sqrt(0.28),8,1800];
	z = flux;

	disp('Number of CT fluxes');
	disp(num2str(size(z)));

	% [x,fval,exitflag,output] = fmincon(@(theta) rml_costfun(theta, z,X,D,E),theta0,[],[],[],[],[0.01 0.01 0.01],[1000 100])

	[x,fval,exitflag,output] = fminsearch(@(theta) rml_costfun(theta, z,X,D,E,ntimes),theta0);

	thetatemp = thetatemp + x; %added 8 times%

	end; % End of timechunk loop timechunk loop
    
    
    
    

	thetatemp = thetatemp ./ 1;  %thereby each month we get one theta, now it's just one%

	disp('Covariance parameter estimate for the month in question');
	disp(thetatemp);
    thetaall = [thetaall; thetatemp];
	end; % End of month or ym loop

	%thetaall = [thetaall; thetatemp ];

	disp('Compiled covariance parameters for all months');
	disp(thetaall)
        save theta_land.mat thetaall
%---------------------------------------------------------------------------------------------------------------%
% END OF SCRIPT
%---------------------------------------------------------------------------------------------------------------%

