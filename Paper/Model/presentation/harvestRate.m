
clf
h = 0:0.05:1;

hold on
axis([-0.1 1.1 0 0.35])
set(gca,'XTickLabelMode','manual')
set(gca,'XTick',[0 0.5 1])
set(gca,'XTickLabel',{'0';'r/2';'r'})

set(gca,'YTickLabelMode','manual')
set(gca,'YTick',[0 0.25 ])
set(gca,'YTickLabel',{'0';'Lh/4'})


plot(h,h.*(1-h),'LineWidth',2.5)
xlabel('h')
ylabel('Harvest Rate')
title('Harvest Rate at Steady State')

print -dpng reducedHarvest.png

clf
hunt = 0.1;

hold on
axis([-0.1 1.1 0 0.35])
set(gca,'XTickLabelMode','manual')
set(gca,'XTick',[0 0.5 0.9 1])
set(gca,'XTickLabel',{'0';'r/2';'r-H';'r'})

set(gca,'YTickLabelMode','manual')
set(gca,'YTick',[0 0.25 ])
set(gca,'YTickLabel',{'0';'Lh/4'})

plot(h,h.*(1-h),'LineWidth',2.5)
plot(h,h.*(1-h-hunt),'r','LineWidth',2.5)
xlabel('h')
ylabel('Harvest Rate')
title('Harvest Rate at Steady State')


print -dpng reducedPopulation.png




