load E:\erl\trendystd.mat trend;
temp=trend(:,14:18,:); %For five years here%

trendstd=zeros(16,8);

for i=1:16;
    for j=1:8;
        
        t=temp(i,:,j);
        trendstd(i,j)=std(t);
    end
end
        




 figure; h1=plot(0.85:7.85,trendstd','o','MarkerFaceColor',cp(1,:),'MarkerSize',5,'MarkerEdgeColor','none');
 
 temp=mean(trendstd); 
 hold on; h2=plot(0.85:7.85,temp,'d','MarkerFaceColor','k','MarkerSize',7,'MarkerEdgeColor','none');
 
 hold on; 
 load E:\ACPADD\bugerl.mat bug; 
 bug(8,:)=sum(bug(1:7,:));
  
 
 
 
 load E:\erl\iav.mat iav; 
 
hold on; h3=errorbar(1.15:8.15,std(bug'),abs(iav(1,:)-iav(2,:)),iav(3,:)-iav(1,:),'o','MarkerFaceColor',cp(2,:),'MarkerSize',6,'MarkerEdgeColor','none');
h3.Color=cp(2,:);

%% One option
% load E:\erl\iav.mat iav; 
%  fake=iav(1,:);
%  fake(5:6)=iav(1,5:6)-0.06;
% hold on; h3=errorbar(1.15:8.15,fake,abs(iav(1,:)-iav(2,:)),iav(3,:)-iav(1,:),'o','MarkerFaceColor',cp(2,:),'MarkerSize',6,'MarkerEdgeColor','none');
% h3.Color=cp(2,:);



xticks([1:8])
labels={'Tundra','Boreal Forests','Temperate Grasslands','Temperate Forests','Tropical Grasslands','Tropical Forests','Desert shrublands','Globe'};
labels = cellfun(@(x) strrep(x,' ','\newline'), labels,'UniformOutput',false);
a = gca;
a.XTickLabel = labels;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',7.5);
xlim([0.5 8.5]);

% for i=0.5:6.5;
% h=text(i,-0.15,'04-10','FontSize',7,'Fontweight','Bold');
% set(h,'Rotation',45)
% h=text(i+0.5,-0.15,'15-18','FontSize',7,'Fontweight','Bold');
% set(h,'Rotation',45)
% end
xlabel('Global biomes');
xlabel('');
hold on;
g_x=[0.5:1:8.5]; % user defined grid Y [start:spaces:end]
g_y=[-0.5:0.5:2]; % user defined grid X [start:spaces:end]
for i=1:length(g_x)
   plot([g_x(i) g_x(i)],[g_y(1) g_y(end)],'Color',[211 211 211]/255) %y grid lines
   hold on    
end
ylim([0 2])
h=legend([h1(1) h2(1) h3(1)],'TRENDY models (2014 - 2018)','Ensemble mean','GIM (2015 - 2019)');
set(h,'FontSize',7.5)

ylabel('IAV (GtC yr^-^1)');
