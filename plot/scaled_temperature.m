%Scaled Temperature Function%

tmin=0; tmax=40; 

topt=20; 
%%topt=15; for tundra and boreal.

t=-5:0.1:45;

ts = (t - tmin).*(t - tmax)./((t - tmin).*(t - tmax) - (t-topt).^2);

ts(t<=0)=0; ts(t>=40)=0;
figure; plot(t,ts,'LineWidth',3);

xlabel('Air Temperature (^oC)');
ylabel('Scaled Temperature');