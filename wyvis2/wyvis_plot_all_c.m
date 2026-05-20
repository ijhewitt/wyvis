function wyvis_plot_all_c(fig,ad,ud)
% plot all channel variables in subplots on large figure
%
% IJH 16 August 2013

     S = [ud.S];
     phi_c = [ud.phi_c];
     Q = [ud.Q];
     M = [ud.M];

    figure(fig); clf; set(gcf,'PaperPositionMode','auto','units','centimeters','Position',[4 4 10 20]);
    set(gcf,'DefaultAxesFontSize',10);
    fontsize = 12;
    linewidth = 1;
    ax1 = axes('position',[.1 .84 .86 .13]);
        plot(ad.x/1e3,ad.phi_b/1e6,'color',0.8*[1 1 1],'linewidth',linewidth); 
        hold on;
        plot(ad.x/1e3,ad.phi_s/1e6,'color',0.8*[1 1 1],'linewidth',linewidth); 
        plot(ad.x/1e3,phi_c/1e6,'b','linewidth',linewidth);
        text(0.5,1,'Hydraulic potential [ MPa ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
        set(gca,'layer','top');
    ax2 = axes('position',[.1 .64 .86 .13]);
        plot(ad.x/1e3,Q,'b','linewidth',linewidth);
        text(0.5,1,'Discharge [ m${}^3$/s ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
    ax3 = axes('position',[.1 .44 .86 .13]);
        plot(ad.x/1e3,S,'b','linewidth',linewidth);
        text(0.5,1,'Cross-section [ m${}^2$ ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
    ax4 = axes('position',[.1 .24 .86 .13]);
        plot(ad.x/1e3,M,'b','linewidth',linewidth);
        xlabel('Distance [ km ]','fontsize',fontsize);
        text(0.5,1,'Melting rate [ m${}^2$/s ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');

end