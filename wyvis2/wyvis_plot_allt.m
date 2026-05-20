function wyvis_plot_allt(fig,td,ud,xi)
% wyvis_plot_allt(fig,td,ud,xi)
% plot all variables over time range of td at pts xi in subplots on large figure
%
% IJH 16 August 2013
     
    if nargin<4, xi = 1:length([ud(1).phi_c]); end
    
     S = [ud.S];
     h = [ud.h];
     phi_c = [ud.phi_c];
     phi = [ud.phi];
     u = [ud.u];
     Q = [ud.Q];
     q = [ud.q];
     N_c = [ud.N_c];
     N = [ud.N];
     M = [ud.M];
     m = [ud.m];
     kappa = [ud.kappa];  
 
    figure(fig); clf; set(gcf,'PaperPositionMode','auto','units','centimeters','Position',[4 4 20 20]);
    set(gcf,'DefaultAxesFontSize',10);
    fontsize = 12;
    linewidth = 1;
    ax1 = axes('position',[.05 .84 .43 .13]);
        plot(td/24/60/60,N_c(xi,:)/1e6,'b','linewidth',linewidth);
        hold on;
        plot(td/24/60/60,N(xi,:)/1e6,'r','linewidth',linewidth);
        text(0.5,1,'Effective pressure [ MPa ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
        set(gca,'layer','top');
    ax2 = axes('position',[.05 .64 .43 .13]);
        plot(td/24/60/60,Q(xi,:),'b','linewidth',linewidth);
        text(0.5,1,'Discharge [ m${}^3$/s ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
    ax3 = axes('position',[.05 .44 .43 .13]);
        plot(td/24/60/60,S(xi,:),'b','linewidth',linewidth);
        text(0.5,1,'Cross-section [ m${}^2$ ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
    ax4 = axes('position',[.05 .24 .43 .13]);
        plot(td/24/60/60,M(xi,:),'b','linewidth',linewidth);
        xlabel('Time [ d ]','fontsize',fontsize);
        text(0.5,1,'Melting rate [ m${}^2$/s ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');

        
    ax6 = axes('position',[.55 .84 .43 .13]);
        plot(td/24/60/60,0*td,'color',0.8*[1 1 1],'linewidth',linewidth); 
        hold on;
        plot(td/24/60/60,kappa(xi,:),'k','linewidth',linewidth); 
        text(0.5,1,'Exchange [ m${}^2$/s ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
        set(gca,'layer','top');
    ax7 = axes('position',[.55 .64 .43 .13]);
        plot(td/24/60/60,q(xi,:),'r','linewidth',linewidth);
        text(0.5,1,'Discharge [ m${}^2$/s ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
    ax8 = axes('position',[.55 .44 .43 .13]);
        plot(td/24/60/60,h(xi,:),'r','linewidth',linewidth);
        text(0.5,1,'Sheet depth [ m ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
    ax9 = axes('position',[.55 .24 .43 .13]);
        plot(td/24/60/60,m(xi,:)*24*60*60*1e3,'r','linewidth',linewidth);
%         xlabel('Time [ d ]','fontsize',fontsize);
        text(0.5,1,'Melting rate [ mm/d ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
    ax10 = axes('position',[.55 .04 .43 .13]);
        plot(td/24/60/60,(365*24*60*60)*u(xi,:),'k','linewidth',linewidth);
        xlabel('Time [ d ]','fontsize',fontsize);
        text(0.5,1,'Ice velocity [ m/y ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');

end

