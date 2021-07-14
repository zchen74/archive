function [ffx]= read_ff(yr);


load /home/grifbake/chen3274/testout_L13/data/landmap.mat landmap;landmap(landmap==2)=0;landmap=reshape(landmap,[],1);
%for yr=2015:2015;
    yer=num2str(yr);

filename=strcat('/home/grifbake/chen3274/odiac/odiac2019_1x1d_',yer,'.nc');end;
ffuel=ncread(filename,'land'); 

temp=zeros(180,360,12);
for i=1:size(ffuel,3),
    temp(:,:,i)=ffuel(:,:,i)';
end
ffuel=temp; clear temp; 


ff=[];

clon=-179.5:179.5;
clat=-79.5:79.5;
clon=-179.5:179.5; clon=repmat(clon,180,1);
clat=-89.5:89.5; clat=repmat(clat',1,360);
load E:\DBAKER\GC_LonLat.mat lon lat; glon=lon; glat=lat;

for ii=1:12;
temp=[];
temp=interp2(clon,clat,ffuel(:,:,ii),glon,glat);
temp=reshape(temp,[],1); temp(landmap~=1)=0;
ff=[ff;temp'];
end



ffx=[];
for i=1:31;ffx=[ffx;ff(1,:)'];end;
for i=1:28;ffx=[ffx;ff(2,:)'];end;
for i=1:31;ffx=[ffx;ff(3,:)'];end;
for i=1:30;ffx=[ffx;ff(4,:)'];end;
for i=1:31;ffx=[ffx;ff(5,:)'];end;
for i=1:30;ffx=[ffx;ff(6,:)'];end;
for i=1:31;ffx=[ffx;ff(7,:)'];end;
for i=1:31;ffx=[ffx;ff(8,:)'];end;
for i=1:30;ffx=[ffx;ff(9,:)'];end;
for i=1:31;ffx=[ffx;ff(10,:)'];end;
for i=1:30;ffx=[ffx;ff(11,:)'];end;
for i=1:31;ffx=[ffx;ff(12,:)'];end; %gc/m2/d%

ffx=ffx./12.*1e6./86400; %umol/m2/s%

save(strcat('ff',yer,'.mat'),'ffx');
clear w p k ff;
%end







