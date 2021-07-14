function [gfed] = read_gfed_hdf5(yr);
%for yr=2015:2015;
    
load /home/grifbake/chen3274/testout_L13/data/landmap.mat landmap;landmap(landmap==2)=0;landmap=reshape(landmap,[],1);
    yer=num2str(yr);
data=[];
for i=1:12,
    
 if yr<=2016;filename=strcat('/home/grifbake/chen3274/GFED4/',yer,'/GFED4.1s_',yer,'.hdf5');end;
 if yr>2016;filename=strcat('/home/grifbake/chen3274/GFED4/',yer,'/GFED4.1s_',yer,'_beta.hdf5');end; 
    if i<10,
    dataname=strcat('/emissions/0',num2str(i),'/C');
    else
       dataname=strcat('/emissions/',num2str(i),'/C');  
    end
    
data(:,:,i)=h5read(filename,dataname); %gC/m2/month%
end

temp=data;
tempp=zeros(720,1440,12);
for i=1:size(temp,3),
    tempp(:,:,i)=temp(:,:,i)';
end

clon=-179.75:0.5:179.75;
clat=-79.75:0.5:79.75;
clon=repmat(clon,360,1);
clat=repmat(clat',1,720);
load E:\DBAKER\GC_LonLat.mat lon lat;
glon=lon; glat=lat;

gf=[];
for i=1:12,
    tmp=[];
tmp=interp2(clon,clat,tempp(:,:,i),lon,lat);
    tmp=reshape(tmp,[],1); tmp(landmap~=1)=0; %we dont care about ocean here%
    gf=[gf;tmp'];
end


gfed=[];
for i=1:31;gfed=[gfed;gf(1,:)'];end;
for i=1:28;gfed=[gfed;gf(2,:)'];end;
for i=1:31;gfed=[gfed;gf(3,:)'];end;
for i=1:30;gfed=[gfed;gf(4,:)'];end;
for i=1:31;gfed=[gfed;gf(5,:)'];end;
for i=1:30;gfed=[gfed;gf(6,:)'];end;
for i=1:31;gfed=[gfed;gf(7,:)'];end;
for i=1:31;gfed=[gfed;gf(8,:)'];end;
for i=1:30;gfed=[gfed;gf(9,:)'];end;
for i=1:31;gfed=[gfed;gf(10,:)'];end;
for i=1:30;gfed=[gfed;gf(11,:)'];end;
for i=1:30;gfed=[gfed;gf(12,:)'];end; %gC/m2/month% 

gfed=gfed/12*1e6/30/86400; %umol/m2/s%
save(strcat('gfed',yer,'.mat'),'gfed');
clear tmp temp tempp;
%end



