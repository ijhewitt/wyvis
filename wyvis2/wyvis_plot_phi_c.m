function wyvis_plot_phi_c(fig,ad,ud)
% plot potential
%
% IJH 16 August 2013
 
    figure(fig); clf; set(gcf,'PaperPositionMode','auto','units','centimeters','Position',[4 4 20 10]);
    set(gcf,'DefaultAxesFontSize',10);
    fontsize = 12;
    linewidth = 2;
    axes('position',[.1 .2 .8 .7]);
        plot(ad.x/1e3,ad.phi_b/1e6,'color',0.8*[1 1 1],'linewidth',linewidth); 
        hold on;
        plot(ad.x/1e3,ad.phi_s/1e6,'color',0.8*[1 1 1],'linewidth',linewidth); 
        plot(ad.x/1e3,[ud.phi_c]/1e6,'b','linewidth',linewidth);
        xlabel('Distance [ km ]','fontsize',fontsize);
        text(0.5,1,'Hydraulic potential [ MPa ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
        set(gca,'layer','top');
    
end