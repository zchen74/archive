%-----------------------------------------------------------------------------------------------------------------------%
% SCRIPT:  inversion_run												%
% PURPOSE: Launch a geostatistical inverse model (GIM), implemented with Lagrange multipliers for non-negativity	%
% S. Miller, Jan. 8, 2016	
% Modified to be used to couple GEOS-Chem forward and adjoint model
% Z. Chen, Dec. 18, 2018
%-----------------------------------------------------------------------------------------------------------------------%

%%
%% please note that here I use Latitude (rows) by Longitude (columns)in a matrix. e.g., 46*72 in the GEOS-Chem 4 (lat) \times 5 (lon) spatial resolution. I Follow this rule and build  E, D, shat0, shat and all matrices possibly related.%

    
    
	% This script will save the estimated fluxes into the following folder
     maxit = 50;
	outpath = './output/';
	
	
	%---------------------------------------------%
	% **Set the covariance matrix parameters**    %
    % Please note the unit for elements in theta  % 
    % are ppm, umol/m2/s, km, and days.           %                             
	%---------------------------------------------%
    
    %NOTE WE ARE NOT USING THE ERROR FOR OBSERVATION MINUS MODEL, AS THE ADJOINT MODEL WILL DO THAT FOR US; I.E., WE ARE NOT USING theta(1) or theta1(1)%
    %but we still leave it here to conform the benchmark code in Github by Dr. Miller%
     theta = [ 1 sqrt(0.31) 1460 5.1 ];  % for land%
     theta1= [ 1 sqrt(0.014) 4678 8.6 ];  %for ocean%
    
	
	% Display the covariance matrix parameters on screen
	disp('Covariance matrix parameters for land');
	disp(theta);
    disp('Covariance matrix parameters for ocean');
	disp(theta1);
	
%---------------------------------%
% **Read in the observations**    %
%---------------------------------%

	disp('Read in the observation vector');

	% Read in the vector of observations here.
	% Subtract off the boundary condition before reading in the observations here.  
	% We'll refer to this vector as "Z".
	% This vector has length (n x 1) where n are the number of observations
	%load ..\data\Z.mat Z
 
%----------------------------%
% **Create the X matrix**    %
%----------------------------%

	disp('Create the X matrix');
	
    %Define n & m%
    %n=length(Z);
    m=3312 * 365; %46 (lat) \times 72 (lon) \times 365 (days)%
    
    
    %load selected environmental drivers X below%
    load ~/testout_2018/shat2015.mat X3; X=X3; clear X3;
    disp('X has been successfully processed'); 

	

%-------------------------%
% Create the E matrix     %
%-------------------------%

	disp('Create the E and D matrices');

	% The E and D matrices are used to construct the Q covariance matrix. 
	% The E matrix describes how covariances in Q decay in space, and the D matrix describes how covariances in Q decay in time.
	% See the following paper by Vineet Yadav for more details on the E and D matrices: http://www.geosci-model-dev.net/6/583/2013/gmd-6-583-2013.html
	% The E matrix has dimensions (r x r)
	% The D matrix has dimensions (q x q)
	

	%----------------------------------------------%
	% **Read in the geographic distance matrix**   %
	%----------------------------------------------%
	
	% This matrix (dimensions r x r) should define the distance between each emissions grid box
	% This matrix is symmetric with zeros on the diagonals. 
	% This matrix should typically has units of km. Whatever units you choose should match the units of theta(3)
	% Refer to Eq. 7 in Gourdji et al., 2010 for detail.
	load /home/grifbake/chen3274/testout_L13/data/deltamat.mat deltamat
	
	
	%------------%
	% Create E   %
	%------------%

	% Spherical covariance model
    disp('Begin to process E with land/ocean indices');
    load /home/grifbake/chen3274/testout_L13/data/landmap.mat; landmap=reshape(landmap,[],1);
    landmap(landmap==2)=nan;%ice is not considered% 
    size(landmap)  
    
    %----------------------------%
% Create the sigmaQ vector   %
%----------------------------%
	
    sigmaQ=zeros(size(landmap),1);
    sigmaQ(landmap==1)=theta(2); %% for land%%
    sigmaQ(landmap==0)=theta1(2); %% for ocean%%
	
    E=zeros(3312,3312);
    for i=1:3312,
        for j=1:3312,
            if landmap(i)+ landmap(j)==2, % when land meets land%
         %       E(i,j)                      = 1 - 1.5 .* (deltamat(i,j) ./theta(3))  + 0.5 .* (deltamat(i,j).^3 ./ theta(3).^3);
         %        temp=E(i,j); temp(deltamat(i,j) > theta(3)) = 0; E(i,j)=temp; %spherical%

               E(i,j) = exp(-1.*(deltamat(i,j)./theta(3)));  %exponential%
              temp=E(i,j); temp(deltamat(i,j) > 3.*theta(3)) = 0; E(i,j)=temp;
              
            elseif landmap(i)+ landmap(j)==0, %when ocean meets ocean%
          %      E(i,j)                      = 1 - 1.5 .* (deltamat(i,j) ./theta1(3))  + 0.5 .* (deltamat(i,j).^3 ./ theta1(3).^3);
          %      temp=E(i,j); temp(deltamat(i,j) > theta1(3)) = 0; E(i,j)=temp;
              
              E(i,j) = exp(-1.*(deltamat(i,j)./theta1(3)));
              temp=E(i,j); temp(deltamat(i,j) > 3.*theta1(3)) = 0; E(i,j)=temp;
            end
        end
    end
    %%note currently we do not consiser the spatial covariance between land and ocean %%  
     clear temp;
      disp('E has been successfully processed'); 
      
    E = (sigmaQ*sigmaQ') .* E;  
      
    % Take the inverse of the E matrix
	Einv  = inv(E);
	
	
%----------------------------%
% Create the D matrix        %
%----------------------------%

	ntimes = m./size(E,1); %% times, i.e, q by q%
	days = 1:ntimes;
	days = days';
	days = days * ones(1,length(days));
	days = abs(days - days');


	%------------%
	% Create D   %
	%------------%
	
    % Option I
	% Spherical covariance model
	% Advantage: it tapers off to zero and is therefore faster to compute with large
	% datasets, e.g., Miller 2018 & 2020;
    %D = 1 - 1.5 .* (days   ./theta(4))  + 0.5 .* (days.^3   ./ theta(4).^3);
    %	D(days   > theta(4)) = 0;
    
    % Option II
    % Exponential, but pay particular attention to the fact that how we
    % define decorrelation length and time  differently between exponential and
    % spherical covariance models. 
    D = exp(-1 .* days ./ theta(4));  D(days > 3.*theta(4)) = 0; Dinv =inv(D);
    D1 = exp(-1 .* days ./ theta1(4));  D1(days > 3.*theta1(4)) = 0; Dinv1 =inv(D1);


%------------------------%
% Create the R matrix    %
%------------------------%
		
% 	disp('Create R');
% 	
% 	% No need to edit this section, unless you want a more complicated setup for R.
% 	R = (theta(1).^2) .* ones(n,1);
% 	R = spdiags(R,0,n,n);


%------------------------------------------------------%
% Create the initial guess for the L-BFGS-B algorithm  %
%------------------------------------------------------%

	disp('Create an initial guess for the L-BFGS-B algorithm');
	
	% The L-BFGS-B algorithm requires an initial guess of the fluxes/emissions. It will then iterate toward the solution.
	% The initial guess can be very important, and a more accurate initial guess is always better than a poor initial guess.
	% The L-BFGS-B will converge quickly on the solution for flux elements that have a large impact on the cost function.
	% However, L-BFGS-B will converge slowly for flux elements that do not have a large impact on the cost function.
	% Here, I will use the deterministic model as the initial guess of the fluxes.

	%--------------------------------------%
	% Fastest but least accurate method    %
	%--------------------------------------%

disp('Load the intial guess of shat')
load ~/testout_2018/shat2015.mat shat0;

count=0;  save /home/grifbake/chen3274/testout_L13/output/count.mat count; clear count;
%-------------------------------------------------%
% Pre-calculate matrix products where possible    %
%-------------------------------------------------%

	disp('Calculate inv(Q) * X');
	% We'll refer to this matrix product as varible "B" from now on
	% B = inv(Q) * X;
    
	B = [];
	m1  = size(E,1);
        landmap=repmat(landmap,1,size(X,2));        

	for j = 1:size(Dinv,1);
	B1 = zeros(m1,size(X,2));
        B2 = zeros(m1,size(X,2));
		for i = 1:size(Dinv,1);
		sel = (m1.*(i-1)+1):(i.*m1)
		B1 = B1 + X(sel,:) .* Dinv(j,i);
                B2 = B2 + X(sel,:) .* Dinv1(j,i);
		end; % End of i loop
        B1(landmap==0)=B2;clear B2;
	temp = Einv * B1; clear B1;
	B = [B; temp];
	end; % End of i loop
	clear B1 temp;

    disp('Calculate several matrix products that will be used in the preconditioner');
	CD = sqrtm(D); CD1=sqrtm(D1);
	CE = sqrtm(E);
	invCD = inv(CD); invCD1=inv(CD1);
	invCE = inv(CE);

        disp('Calculate inv(sqrtm(Q)) * X');   
        A = [];

        for j = 1:ntimes;
        A1 = zeros(m1,size(X,2));
        AA1=zeros(m1,size(X,2));
                for i = 1:size(Dinv,1);
                sel = (m1.*(i-1)+1):(i.*m1);
                A1 = A1 + X(sel,:) .* invCD(j,i);
                AA1 = AA1 + X(sel,:) .* invCD1(j,i); 
                end; % End of i loop
        A1(landmap==0)=AA1; clear AA1;  
        temp = invCE * A1;  clear A1;
        A = [A; temp];
        end; % End of i loop
        clear A1 temp;	
   

        landmap=landmap(:,1);
 %transform shat0 to shat0*%
        Qx=[];
        for j = 1:ntimes;
        Qx1   = zeros(m1,1);
        QQx1  =zeros(m1,1);
               for i = 1:ntimes;
               sel = (m1.*(i-1)+1):(i.*m1);
                Qx1 = Qx1 + shat0(sel) .* invCD(j,i);
                QQx1 = QQx1 + shat0(sel) .* invCD1(j,i);
                end; % End of i loop
        Qx1(landmap==0)=QQx1; clear QQx1;
        temp = invCE * Qx1;   clear QX1;
        Qx   =  [Qx; temp];
        end; % End of i loop
        clear Qx1 temp;
        shat0=Qx; clear Qx;
        disp('size of initial guess'); 
        size(shat0)	
%--------------------------------------------------%
% Estimate the fluxes using Lagrange multipliers   %
%--------------------------------------------------%

	disp('Run inversion (optimized with the L-BFGS algorithm)');
%-------------------------------------------------------------%
% Create function handles for the cost function and gradient  %

	f1 = @(shat) cost_gradient_fun(X, A, B, Einv, Dinv, Dinv1, CD, CD1, CE, shat, landmap); 

        % Create an empty flux estimate
        shat = [];

        options = struct('HessUpdate','lbfgs','GradObj','on','Display','iter','MaxIter',maxit,'GradConstr',false);

        [shat,costfun,exitflag,gradient] = fminlbfgs(f1,shat0,options);

%-----------------------------------------------%
% Transform fluxes back to normal space         %
%-----------------------------------------------%

	stemp = [];
 
 	ntimes = size(D,1);
         for j = 1:ntimes;
         Qx1   = zeros(m1,1);
         QQx1  = zeros(m1,1);
                 for i = 1:ntimes;
                 sel = (m1.*(i-1)+1):(i.*m1);
                Qx1 = Qx1 + shat(sel) .* CD(j,i);
               QQx1 = QQx1 + shat(sel) .* CD1(j,i);
                end; % End of i loop
        Qx1(landmap==0)=QQx1; clear QQx1;      
        temp =  CE * Qx1; clear Qx1;
        stemp   =  [stemp; temp];
        end; % End of i loop
        clear Qx1 temp;

	shat = stemp;
	
	
%------------------------------%
% Write the outputs to file    %
%------------------------------%

         disp('Writing outputs to file');
        outname = strcat(outpath,'fluxes_LBFGS',num2str(maxit),'.csv');
        dlmwrite(outname,full(shat),',');

        disp('Outputs written to file');
        disp(outname);
	


%-----------------------------------------------------------------------------------------------------------------------%
% END OF SCRIPT														%
%-----------------------------------------------------------------------------------------------------------------------%
