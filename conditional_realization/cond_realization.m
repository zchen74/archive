

%load the budget%
load E:\ACPADD\bugerl.mat bug;
budget4=bug';
budget4(:,8)=sum(budget4(:,1:7),2);

% budget4 will have fivei (years) rows and eight columns (seven biomes plus the globe) 

%load posterior errors for the seven biomes and the globe.
err=reshape(err,1,[]);


iav=zeros(6,8); 

for i=1:8;
    i
 tt=[];   
for j=1:100000; %likely takes 10 sec%
temp=zeros(4,1);
temp(1,1)=randn*er(1,i)./1.96+budget4(1,i); %posterior errors are w/ 95% confidence, so we divide by 1.96%
temp(2,1)=randn*er(2,i)./1.96+budget4(2,i);
temp(3,1)=randn*er(3,i)./1.96+budget4(3,i);
temp(4,1)=randn*er(4,i)./1.96+budget4(4,i);
temp(5,1)=randn*er(5,i)./1.96+budget4(5,i);
tt=[tt;max(temp)-min(temp)];%get the range -- max -min%
end

iav(1,i)=quantile(tt,0.17);
iav(2,i)=quantile(tt,0.83);
iav(3,i)=quantile(tt,0.05);
iav(4,i)=quantile(tt,0.95);
iav(5,i)=quantile(tt,0.025);
iav(6,i)=quantile(tt,0.975);
end
iav
