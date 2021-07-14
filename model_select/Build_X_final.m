for yr=2015:2018;
x=read_merra2(yr);
gfed=read_gfed_hdf5(yr);
ffx=read_ff(yr);
load ocn.mat ocn; %from NASA ECCO DARWIN; and alreadt filted by biome%
yer=num2str(yr);
disp(yer);
X=[];



load /home/grifbake/chen3274/testout_L13/data/landmap.mat landmap;landmap(landmap==2)=0;landmap=reshape(landmap,[],1);
landmap=repmat(landmap,1,364);landmap=reshape(landmap,[],1);
load /home/grifbake/chen3274/testout_L13/data/biome_scot.mat biome;biome=reshape(biome,[],1);
biome=repmat(biome,1,364);biome=reshape(biome,[],1);
load /home/grifbake/chen3274/testout_L13/data/GC_LonLat.mat lon lat; LAT=repmat(lat,1,1,364); LAT=reshape(LAT,[],1);

%Intercepts for land%
for i=1:7;
    temp=zeros(3312*364,1);
    temp(biome==i)=1;
    X=cat(2,X,temp);
end

%Intercepts for ocean%
for i=1;
    temp=zeros(3312*364,1);
    temp(biome==0)=1;
    X=cat(2,X,temp);
end

num= xxxx; %how many drivers are there that you include%
XX=[]; %now we have 7x clms%
for i=1:num;
    for j=1:7;
        temp=x(:,i);
        temp(biome~=j)=nan;
        temp=(temp-nanmean(temp))./nanstd(temp);
        temp(isnan(temp))=0;
        XX=cat(2,XX,temp);
    end
end
X=cat(2,X,XX);

%scaled temp%
for i=3;
    for j=1:7;
        temp=x(:,i);
        temp(biome~=j)=nan;
        
     if j<=2; opt=15; end;
     if j>2; opt=20; end;
        
        temp = temp.*(temp-40)./(temp.*(temp-40)-(temp-opt).^2);
        temp=(temp-nanmean(temp))./nanstd(temp);
        temp(isnan(temp))=0;
        X=cat(2,X,temp);
    end
end

X=cat(2,X,gfed+ffx+ocn);

save(strcat('/panfs/roc/groups/8/grifbake/chen3274/testout_',yer,'/X.mat'),'X');
end
