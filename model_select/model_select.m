%Model selection on candiate environmental drivers; i.e., intercepts, env dirvers, fossil fuel+biomass burning+ocean%


disp('load HX and y for each individual year')
yr=2015;

load(strcat('HX',num2str(yr),'.mat'),'HX','y'); 

X=HX;

p=size(HX,2)-7; %7 intercepts% 
theta=1.193;
n=size(y,1);
nstar=1e4;
n=double(n); nstar=double(nstar); 

A = y' * y; A=(1./theta).*A;
B = y' * X; B=(1./theta).*B;
C = X' * X; C=(1./theta).*C;

detpsi=0; 

%save write_input_for_selection.mat y A B C p detpsi nstar

[ selmodel beta1 ] = heuristic_bb_zc(y, A, B, C, p, detpsi, nstar);
disp('The number of selected indice are')
disp(num2str(size(selmodel,2)));


%The idea is that we run the model selection for each year, and use the drivers that are always selected%

