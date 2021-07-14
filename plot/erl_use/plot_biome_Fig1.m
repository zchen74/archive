
load biome.mat biome;
load GC_LonLat.mat; LON=lon(1,:)';LAT=lat(:,1);
%load landmap.mat; landmap(landmap==2)=0; 
%LON=-177.5:5:177.5; LAT=-88:4:88;
%figure;
data=biome;
data(biome==0)=nan;
data=8-data;
imagesc(LON,LAT,data,'AlphaData',~isnan(data));set(gca,'YDir','normal');colorbar; 
%contourf(LON',LAT,data(1:45,:));colorbar; 
%caxis([1 7]);
%caxis([-1 1])
%colormap(brewermap(7,'Paired'));
 %set(colorbar,'XTick',[1:1:7]); 
%%
% [hC hC]=contourf(LON,LAT,data);colorbar;
% set(hC,'LineStyle','none');
colormap(brewermap(7,'Paired'));
set(colorbar,'XTick',[1:1:7]); 
caxis([1 7]);
%hold on; load coast; plot(long,lat,'Color',[17 17 17]/255);

hold on;
geoshow('landareas.shp', 'FaceColor', 'none');
xlabel('LON'); ylabel('LAT');% title(num2str(i)
   
%legend({'Tundra','Boreal Forests','Temperate Grasslands','Temperate Forests','Tropical Grasslands','Tropical Forests','Desert and shrublands'},'FontSize',8);

%ocean_mask_function



c=colorbar('southoutside');
c.TickLabels={'Tundra','Boreal Forests','Temperate Grasslands','Temperate Forests','Tropical Grasslands','Tropical Forests','Desert shrublands'};
temp=c.TickLabels;
colorbar off; 
temp1=temp;

for i=1:7;
    temp1(i)=temp(8-i);
end
colorbar off;
hold on;
c=colorbar('southoutside');
%labels = {'line1 line2','line1 line2','line1 line2'};
c.TickLabels = cellfun(@(x) strrep(x,' ','\newline'), c.TickLabels,'UniformOutput',false);
c.TickLabels = cellfun(@(x) strrep(x,' ','\newline'), temp1,'UniformOutput',false);
c.FontSize=8;
set(gca,'position',[0.11 0.21 0.7750 0.75])
set(c,'position',[ 0.11    0.0729    0.7750    0.0308])
caxis([0.5 7.5]);



xlabel('LON','FontSize',10,'FontWeight','Bold');
ylabel('LAT','FontSize',10,'FontWeight','Bold');
% text(-179,86,'(b)')
% saveas(gcf,'Olson.epsc');
% saveas(gcf,'Olson.pdf');
% saveas(gcf,'Olson.tif')
