load /home/grifbake/chen3274/testout_L13/data/landmap.mat landmap;
     landmap(landmap==2)=0; %0 for ocean and 1 for land, 2 is for ocean%
    
    lmp=reshape(landmap,[],1); %reshape for specific purpose that I cannot tell you%
     % Create the spatial distance matrix
    load('/home/grifbake/chen3274/testout_L12/data/deltamat.mat','deltamat');




count=[];
for i=1:3312,
    if lmp(i)==0,  %all you need to change%
     count=[count;i];   
    end
end
distmat=zeros(length(count),length(count));

for i=1:length(count),
    for j=1:length(count),
        distmat(i,j)=deltamat(count(i),count(j));
    end
end


save distmat_ocean.mat distmat 


