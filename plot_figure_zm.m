function  plot_figure_zm(lat,lev,var,mymap,coeff,xlim_value)
figure
contourf(lat,lev,var',-coeff*20:coeff/10:coeff*20,'Linecolor','none')
set(gca, 'YDir','reverse')
xlabel('Latitude')
ylabel('Pressure (hPa)')
set(gca,'YScale','log');
xlim(xlim_value)
%xlim([-19 61])
%xlim([-80 -20])
ylim([50 1000])
%ylim([100 600])
%ylim([10 850])
xticks(-90:10:90)
%xticklabels({'90S','75S','60S','45S','30S','15S','EQ','15N','30N','45N','60N','75N','90N'})
xticklabels({'90S','80S','70S','60S','50S','40S','30S','20S','10S',...
    'EQ','10N','20N','30N','40N','50N','60N','70N','80N','90N'})
yticks([10 20 30 50 70 100 150 200 250 300 400 500 700 850 1000])
yticklabels({'10','20','30','50','70','100','150','200','250','300','400','500','700','850','1000'})
clim([-coeff coeff]);
colormap(mymap)
set(gca,'fontsize',8)
line([0, 0], ylim, 'Color', [0.8 0 0], 'LineWidth', 1);