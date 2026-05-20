function wyvis_plot_all(fig,ad,ud)
% wyvis_plot_all(fig,ad,ud)
% plot all variables in subplots on large figure
%
% IJH 16 August 2013

     S = [ud.S];
     h = [ud.h];
     phi_c = [ud.phi_c];
     phi = [ud.phi];
     u = [ud.u];
     Q = [ud.Q];
     q = [ud.q];
     M = [ud.M];
     m = [ud.m];
     kappa = [ud.kappa];  
     Sw = [ud.Sw];
     hw = [ud.hw];
 
    figure(fig); clf; set(gcf,'PaperPositionMode','auto','units','centimeters','Position',[4 4 20 20]);
    set(gcf,'DefaultAxesFontSize',10);
    fontsize = 12;
    linewidth = 1;
    ax1 = axes('position',[.05 .84 .43 .13]);
        plot(ad.x/1e3,ad.phi_b/1e6,'color',0.8*[1 1 1],'linewidth',linewidth); 
        hold on;
        plot(ad.x/1e3,ad.phi_s/1e6,'color',0.8*[1 1 1],'linewidth',linewidth); 
        plot(ad.x/1e3,phi_c/1e6,'b','linewidth',linewidth);
        plot(ad.x/1e3,phi/1e6,'r','linewidth',linewidth);
        text(0.5,1,'Hydraulic potential [ MPa ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
        set(gca,'layer','top');
    ax2 = axes('position',[.05 .64 .43 .13]);
        plot(ad.x/1e3,Q,'b','linewidth',linewidth);
        hold on; plot(ad.x/1e3,Q+ad.W.*q,'k','linewidth',linewidth); % total discharge
        text(0.5,1,'Discharge [ m${}^3$/s ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
    ax3 = axes('position',[.05 .44 .43 .13]);
        plot(ad.x/1e3,S,'b','linewidth',linewidth);
%         hold on; plot(ad.x/1e3,Sw,'b','linewidth',linewidth); 
        text(0.5,1,'Cross-section [ m${}^2$ ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
    ax4 = axes('position',[.05 .24 .43 .13]);
        plot(ad.x/1e3,M,'b','linewidth',linewidth);
%         xlabel('Distance [ km ]','fontsize',fontsize);
        text(0.5,1,'Melting rate [ m${}^2$/s ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
    ax5 = axes('position',[.05 .04 .43 .13]);
        plot(ad.x/1e3,ad.Z_b,'k','linewidth',linewidth); 
        hold on;
        plot(ad.x/1e3,ad.Z_s,'k','linewidth',linewidth); 
        xlabel('Distance [ km ]','fontsize',fontsize);
        text(0.5,1,'Elevation [m]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');

    ax6 = axes('position',[.55 .84 .43 .13]);
        plot(ad.x/1e3,0*ad.x,'color',0.8*[1 1 1],'linewidth',linewidth); 
        hold on;
        plot(ad.x/1e3,kappa,'k','linewidth',linewidth); 
        text(0.5,1,'Exchange [ m${}^2$/s ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
        set(gca,'layer','top');
    ax7 = axes('position',[.55 .64 .43 .13]);
        plot(ad.x/1e3,q,'r','linewidth',linewidth);
        text(0.5,1,'Discharge [ m${}^2$/s ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
    ax8 = axes('position',[.55 .44 .43 .13]);
        plot(ad.x/1e3,h,'r','linewidth',linewidth);
%         hold on; plot(ad.x/1e3,hw,'r','linewidth',linewidth); 
        text(0.5,1,'Sheet depth [ m ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
    ax9 = axes('position',[.55 .24 .43 .13]);
        plot(ad.x/1e3,m*24*60*60*1e3,'r','linewidth',linewidth);
%         xlabel('Distance [ km ]','fontsize',fontsize);
        text(0.5,1,'Melting rate [ mm/d ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
    ax10 = axes('position',[.55 .04 .43 .13]);
        plot(ad.x/1e3,u*365*24*60*60,'k','linewidth',linewidth);
        xlabel('Distance [ km ]','fontsize',fontsize);
        text(0.5,1,'Ice velocity [ m/y ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');

end