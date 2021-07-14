load conttrend.mat cont

%bar(1:17,cont(1:16,2)*100,'b');
%hold on; bar(16,cont(16,2)*100,'r');
[temp,I]=sort(cont(1:17,2)*100);
labels={'CABLE-POP','CLASS-CTEM','CLM5.0','DLEM','JSBACH','JULES','LPJ','LPJ-GUESS','LPX','OCN','ORCHIDEE','ORCHIDEE-CNP','SDGVM','SURFEX','VISIT','ISAM','This study'};
labels2=labels;
labels2=labels(I);
figure; 
bar(1:17,temp,'FaceColor',cp(1,:));
%text(11,83,'Tropical grasslands','Fontweight','bold');
ylabel('Contribution (%)')
xticks([1:17])
a = gca;
a.XTickLabel = labels2;
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',7.5)
xlim([0.4 17.6])
xtickangle(45)
xlabel('Models')
hold on; bar(14,temp(14),'FaceColor',cp(2,:));
xlabel('Models')
xlabel('');
ylim([-10 75])
%text(0.5,83,'(a)')
title('(a) Tropical grasslands','FontSize',9)


load conttrend.mat cont

%bar(1:17,cont(1:16,2)*100,'b');
%hold on; bar(16,cont(16,2)*100,'r');

[temp,I]=sort(cont(1:17,3)*100);
labels={'CABLE-POP','CLASS-CTEM','CLM5.0','DLEM','JSBACH','JULES','LPJ','LPJ-GUESS','LPX','OCN','ORCHIDEE','ORCHIDEE-CNP','SDGVM','SURFEX','VISIT','ISAM','This study'};
%labels={'BIOME-BGC','CLASS-CTEM','CLM4VIC','CLM','DLEM','GTEC','ISAM','LPJ','ORCHIDEE','SIB3','SIBCASA','TEM6','TRIPLEX','VEGAS','VISIT','This study'};
labels2=labels;
labels2=labels(I);
figure; 
bar(1:17,temp,'FaceColor',cp(1,:));
%text(11,83,'Tropical forests','Fontweight','bold');
ylabel('Contribution (%)')
xticks([1:17])
a = gca;
a.XTickLabel = labels2;
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',7.5)
xlim([0.4 17.6])
xtickangle(45)
xlabel('Models')
hold on; bar(9,temp(9),'FaceColor',cp(2,:));
xlabel('Models')
xlabel('');
ylim([-10 75])
%text(0.5,83,'(b)')
title('(b) Tropical forests','FontSize',9)



load conttrend.mat cont

%bar(1:17,cont(1:16,2)*100,'b');
%hold on; bar(16,cont(16,2)*100,'r');
[temp,I]=sort(cont(1:17,1)*100);
labels={'CABLE-POP','CLASS-CTEM','CLM5.0','DLEM','JSBACH','JULES','LPJ','LPJ-GUESS','LPX','OCN','ORCHIDEE','ORCHIDEE-CNP','SDGVM','SURFEX','VISIT','ISAM','This study'};
%labels={'BIOME-BGC','CLASS-CTEM','CLM4VIC','CLM','DLEM','GTEC','ISAM','LPJ','ORCHIDEE','SIB3','SIBCASA','TEM6','TRIPLEX','VEGAS','VISIT','This study'};
labels2=labels;
labels2=labels(I);
figure; 
bar(1:17,temp,'FaceColor',cp(1,:));
%text(12,83,'Extra-tropics');
ylabel('Contribution (%)')
xticks([1:17])
a = gca;
a.XTickLabel = labels2;
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',7.5)
xlim([0.4 17.6])
xtickangle(45)
xlabel('Models')
hold on; bar(5,temp(5),'FaceColor',cp(2,:));
xlabel('Models')
ylim([-10 75])
%text(0.5,83,'(c)')
title('(c) Extra-tropics','FontSize',9)

