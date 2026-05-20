function wyvis_plot_phit(fig,td,ud,xi)
% plot ice velocity over time range of td at pts xi
%
% IJH 9 Sept 2016

    if nargin<4, xi = 1:length([ud(1).u]); end

    u = [ud.u];
  
    figure(fig); clf; set(gcf,'PaperPositionMode','auto','units','centimeters','Position',[4 4 20 20]);
    set(gcf,'DefaultAxesFontSize',10);
    fontsize = 12;
    linewidth = 2;
    axes('position',[.1 .2 .8 .7]);
        plot(td/24/60/60,(365*24*60*60)*[ud.u],'k','linewidth',linewidth);
        xlabel('Time [ d ]','fontsize',fontsize);
        text(0.5,1,'Ice velocity [ m/y ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
        set(gca,'layer','top');
    
end

