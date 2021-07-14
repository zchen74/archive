%We use scatter_kde to show the color-coded plot%
%Scatter plot where each point is colored by density%
figure; scatter_kde(y3, y2, 'filled', 'MarkerSize', 10);

ylim([385 420]);xlim([385 420]);
hold on; plot(385:0.1:420,385:0.1:420,'k')
xlabel('observed XCO_2 (ppm)');
ylabel('modeled XCO_2 (ppm)')
text(386,419,'bias = -0.08 ppm','Fontsize',12);
text(386,416.5,'RMSE = 1.10 ppm','Fontsize',12);
title('(a) Year 2015','Fontsize',14);
set(gca,'xtick',[]);
