%-----------------------------------------------------------------------------------------------------------------------------------------------%
% FUNCTION; cost_gradient_fun
% testout_L13
% PURPOSE:  Compute the geostatistical inversion cost function. Use the Kronecker product to make these calculations more effeicient.		%
%	     This script will calculate both the cost function and the gradient (when the latter is requested).					%
% S. Miller, Sept. 21, 2018																%
% Modified for the GEOS-Chem adjoint on-line modeling runnning
% Z.Chen, Jan 1st, 2019
%-----------------------------------------------------------------------------------------------------------------------------------------------%


%-------------%
% NOTES: Make sure in each idl script, put the 'exit' command to close the idl interface%
%-------------%


function [ f,g ] = cost_gradient_fun( X, A, B, Einv, Dinv, Dinv1, CD, CD1, CE, shat, landmap);

load /home/grifbake/chen3274/testout_L13/output/count.mat count;
count=count+1; %here count does not necesarily equal to the iterations, and it does not have to%
save /home/grifbake/chen3274/testout_L13/output/count.mat count;

%save the estimated shat at each count%
save(strcat('/home/grifbake/chen3274/testout_L13/output/shat',num2str(count),'.mat'),'shat'); 

% Important inputs
ntimes = size(Dinv,1);
m =3312 * 365;
m1 = m ./ ntimes; % Divide m by the total number of time periods in the inversion


%---------------------------------------%
%***** CALCULATE THE COST FUNCTION *****%
%---------------------------------------%

disp('Calculating the cost function');
tic;

%%
%% Step 1. Write emissions (i.e., Matlab variable "shat") to a binary punch file
%	GEOS-Chem will then read in that bunch file and use those emissions.
%	Note: This section has been debugged and works.

% Calculate chol(Q)*shat
% I've included a transpose here.
CDt = CD';
CDt1=CD1';
Qx = [];

for j = 1:ntimes;
    Qx1   = zeros(m1,1);
    QQx1  = zeros(m1,1);
    for i = 1:ntimes;
        sel = (m1.*(i-1)+1):(i.*m1);
        Qx1 = Qx1 + shat(sel) .* CDt(j,i);
        QQx1= QQx1 + shat(sel).* CDt1(j,i);
    end; % End of i loop
    Qx1(landmap==0)=QQx1; clear QQx1;
    temp =  CE' * Qx1; clear Qx1;
    Qx   =  [Qx; temp];
end; % End of i loop
clear Qx1 temp;

disp('Now we estimate beta1');
B2 = [];
for j = 1:ntimes;
    A1 = zeros(m1,1);
    AA1 = zeros(m1,1);
    for i = 1:ntimes;
        sel = (m1.*(i-1)+1):(i.*m1);
        A1 = A1 + Qx(sel) .* Dinv(j,i);
        AA1= AA1 + Qx(sel).* Dinv1(j,i);
    end; % End of i loop
A1(landmap==0)=AA1; clear AA1;
    temp = Einv * A1;
    B2 = [B2; temp];
end; % End of i loop
clear A1 temp;

% Estimate the drift coefficients using the current estimate for the fluxes
beta1 = (X' * B) \ (X' * B2);
clear B2;
disp('Estimated drift coefficients');
disp(num2str(beta1));
save(strcat('/home/grifbake/chen3274/testout_L13/output/bta',num2str(count),'.mat'),'beta1');


%To pass X through h(), we need X to be in bpch format, as described below%
%we convert X to netcdf format first, and then convert netcdf to bpch format%
%% WHY DO NOT I DIRECTLY CONVERT X TO BPCH FORMAT--BECAUSE I REALLY DONT KNOW HOW%%

temp=Qx.*6.02e+13; %convert to the unit of molec/cm2/s, the unit needed from GC model%
temp=reshape(temp,46,72,ntimes); %again here, 46(lat) \times 72(lon) \times days%
for i=1:ntimes,
    ncid = netcdf.open(strcat('/panfs/roc/groups/8/grifbake/chen3274/ntei/zc46/zichong',num2str(i),'.nc'),'WRITE'); %where I save the netcdf-format emission input%
    varid = netcdf.inqVarID(ncid,'CO2SRCE');
    netcdf.putVar(ncid,varid,temp(:,:,i)'); 
    netcdf.close(ncid);
end
clear temp; 

%now let us call the use of IDL software, via the unix command%
%convert the netcdf format to bpch format%
unix('/soft/idl/8.0/bin/idl /home/grifbake/chen3274/ntei/zc46/run_final.pro'); 
%%

%% Step 2. run the GC forward and GC adjoint model, to get Hs and the gradient

%submit the model run -- forward and adjoint%
unix('sbatch /home/grifbake/chen3274/bbyrne/runs/v8-02-01/geos9/run.4')


%below is the command to monitor if the GC model is still runing or has completed; it works for Minnesota supercomputer, and please adjust it to faciliate MARCC%

%The only purpose is to check if the model is still running or completed%
[status,cmdout]=unix('sq'); %check the status, please adjust to whatever alias that you use in your supercomputer%
while cmdout(length(cmdout)-17)~='C'; %C means complete%
    disp('GEOS-Chem forward and adjoint model is still running');
    pause(180); %pause for three minutes%
    [status,cmdout]=unix('sq');
end; % 
disp('GEOS-Chem model finishes running, now start data processing')

% Convert the GEOS-Chem outputs from binary punch file to netcdf
unix('/soft/idl/8.0/bin/idl /home/grifbake/chen3274/bbyrne/runs/v8-02-01/geos9/OptData/run.pro')



filename='/home/grifbake/chen3274/bbyrne/runs/v8-02-01/geos9/OptData/gdt.csv';
adj=csvread(filename); adj=reshape(adj,72,46,[]); adj=adj.*6.02.*(10^13); %(1/(umol/m2/s))%
adjj=[];
for i=1:size(adj,3),
    adjj(:,:,i)=adj(:,:,i)'; %transpose the matrix%
end

adj=reshape(adjj,[],1);
clear adjj;

temp=dlmread('/home/grifbake/chen3274/bbyrne/runs/v8-02-01/geos9/OptData/cfn.01');
L1=temp(2); clear temp;
disp('Cost function value of L1');
disp(num2str(L1));

%------------------------------%

Gshat = shat - A * ((X' * B) \ (A' * shat));
L2 = shat' * Gshat;
clear Gshat;

disp('Cost function value of L2');
disp(num2str(L2));


%-------------------------------------------------%
% Add up the two components of the cost function  %
%-------------------------------------------------%

f = L1+ L2;
%f = L1;
disp('Time elapsed for cost function calculations');
disp(toc);

disp('Cost function value');
disp(f);

save(strcat('/home/grifbake/chen3274/testout_L13/output/cfun',num2str(count),'.mat'),'L1','L2','f');
%---------------------------------------%
%***** CALCULATE THE GRADIENT      *****%
%---------------------------------------%

if ( nargout > 1 );
    
    disp('Calculating the gradient');
    
    %------------------------------------------------------%
    % Calculate the observation component of the gradient  %
    %------------------------------------------------------%
    L1 = -adj; %to rogorously match how Scot builds it, note in lines below how we write g%
    

   %%%% L1 = -adj; %now we have to test%
    % Multiply the variable "L1" by chol(Q)
    Qx = [];
    
    for j = 1:ntimes;
        Qx1   = zeros(m1,1);
        for i = 1:ntimes;
            sel = (m1.*(i-1)+1):(i.*m1);
            Qx1 = Qx1 + L1(sel,:) .* CD(j,i);
        end; % End of i loop
        temp =  CE * Qx1;
        Qx   =  [Qx; temp];
    end; % End of i loop
    clear Qx1 temp;
    L1 = Qx;
    
    
    %------------------------------------------------%
    % Calculate the prior component of the gradient  %
    %------------------------------------------------%
    
    % Calculate G * shat
    % L2 = [I + (sQinv * X) * ((X' * Qinv * X) \ (sQinv * X)') ] * shat
    
    L2 = shat - A * ((X' * B) \ (A' * shat));
    
    %--------------------------------------------%
    % Add up the two components of the gradient  %
    %--------------------------------------------%
    
    g = L2 - L1;
    %note L1 I get from the GC adjoint model is negative%
    disp('Time elapsed for gradient calculations');
    disp(toc);
    
    disp('Gradient mean of L1');
    disp(num2str(mean(L1)));
    
    disp('Gradient mean of L2');
    disp(num2str(mean(L2)));
    
    
    
    
end; % End of nargout if statement

%save(strcat('/home/grifbake/chen3274/testout_L10/output/shat',num2str(count),'.mat'),'shat');
%-----------------------------------------------------------------------------------------------------------------------------------------------%
% END OF FUNCTION
%-----------------------------------------------------------------------------------------------------------------------------------------------%
