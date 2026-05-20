function wyvis_plot_phi_ct(fig,td,ud,xi)
% plot potential over time range of td at pts xi
%
% IJH 16 August 2013

    if nargin<4, xi = 1:length([ud(1).phi_c]); end

    phi_c = [ud.phi_c];
  
    figure(fig); clf; set(gcf,'PaperPositionMode','auto','units','centimeters','Position',[4 4 20 20]);
    set(gcf,'DefaultAxesFontSize',10);
    fontsize = 12;
    linewidth = 2;
    axes('position',[.1 .2 .8 .7]);
        plot(td/24/60/60,phi_c(xi,:)/1e6,'b','linewidth',linewidth);
        xlabel('Time [ d ]','fontsize',fontsize);
        text(0.5,1,'Hydraulic potential [ MPa ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
        set(gca,'layer','top');
    
end

