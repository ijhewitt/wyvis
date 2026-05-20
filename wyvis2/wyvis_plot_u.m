function wyvis_plot_phi(fig,ad,ud)
% plot ice velocity
%
% IJH 9 Sept 2016
 
    figure(fig); clf; set(gcf,'PaperPositionMode','auto','units','centimeters','Position',[4 4 20 10]);
    set(gcf,'DefaultAxesFontSize',10);
    fontsize = 12;
    linewidth = 2;
    axes('position',[.1 .2 .8 .7]);
        plot(ad.x/1e3,(365*24*60*60)*[ud.u],'k','linewidth',linewidth);
        xlabel('Distance [ km ]','fontsize',fontsize);
        text(0.5,1,'Ice velocity [ m/y ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
        set(gca,'layer','top');
    
end