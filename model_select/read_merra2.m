function [x] = read_merra2(yr);

load /home/grifbake/chen3274/testout_L13/data/landmap.mat landmap;
landmap(landmap==2)=0; %0 for ocean and 1 for land, 2 is for ocean%
landmap=reshape(landmap,[],1);
%For year 2015%
%Written by Zichong Chen, Johns Hopkins University%
%Feb 22nd, 2019%
%for yr=2015;

    yer=num2str(yr);
mm=[31 28 31 30 31 30 31 31 30 31 30 30];

prec=[]; airt=[]; lai=[]; sw=[];par=[];
ts=[];sm=[];sh=[]; %soil temp,soil moisture, specific humidity %
for i=1:12,
  for j=1:mm(i),

      if i<10,
          ii=strcat('0',num2str(i));
      else ii=num2str(i);
      end
      if j<10,
          jj=strcat('0',num2str(j));
      else jj=num2str(j);
      end

filename=strcat('/home/grifbake/chen3274/shared/data/geos-chem/ExtData/GEOS_4x5/MERRA2/',yer,'/',ii,'/MERRA2.',yer,ii,jj,'.A1.4x5.nc4');
%filename='MERRA2.20150701.A1.4x5.nc4';
temp1=ncread(filename,'PRECTOT'); 
temp1=nanmean(temp1,3).*3600.*24; %mm/day%/;
temp1=reshape(temp1',[],1);  %follow my preference here%
temp1(landmap==0)=0; %we dont care about ocean%

temp2=ncread(filename,'T2M'); 
temp2=nanmean(temp2,3)-273.15; %oC%
temp2=reshape(temp2',[],1);  %follow my preference here%
temp2(landmap==0)=0; %we dont care about ocean%

temp3=ncread(filename,'PARDF')+ncread(filename,'PARDR'); 
temp3=nanmean(temp3,3); %unitless%;
temp3=reshape(temp3',[],1);  %follow my preference here%
temp3(landmap==0)=0; %we dont care about ocean%

temp4=ncread(filename,'SWGDN'); 
temp4=nanmean(temp4,3); %W/m2%;
temp4=reshape(temp4',[],1);  %follow my preference here%
temp4(landmap==0)=0; %we dont care about ocean%

temp5=ncread(filename,'TS');%soil temp%
temp5=nanmean(temp5,3)-273.15; %oC%
temp5=reshape(temp5',[],1);  %follow my preference here%
temp5(landmap==0)=0; %we dont care about ocean%

temp6=ncread(filename,'GWETTOP'); %soil wetness%
temp6=nanmean(temp6,3); %v/v%
temp6=reshape(temp6',[],1);  %follow my preference here%
temp6(landmap==0)=0; %we dont care about ocean%

temp7=ncread(filename,'QV2M'); %specific humidity%
temp7=nanmean(temp7,3); %V/V%
temp7=reshape(temp7',[],1);  %follow my preference here%
temp7(landmap==0)=0; %we dont care about ocean%

prec=[prec; temp1]; airt=[airt; temp2]; par=[par;temp3]; sw=[sw;temp4];ts=[ts;temp5];sm=[sm;temp6];
sh=[sh;temp7];
clear temp1 temp2 temp3 temp4 temp5 temp6 temp7;

  end
end

x=cat(2,prec,airt,par,sw,ts,sm,sh);
%There needs some work to convert to and get scaled temperature, monthly averaged precipitation, and relative humidity%


