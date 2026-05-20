function wyvis_plot_all_ct(fig,td,ud,xi)
% plot all channel variables over time range of td at pts xi in subplots on large figure
%
% IJH 16 August 2013
     
    if nargin<4, xi = 1:length([ud(1).phi_c]); end
    
     S = [ud.S];
     phi_c = [ud.phi_c];
     Q = [ud.Q];
     N_c = [ud.N_c];
     M = [ud.M];
 
    figure(fig); clf; set(gcf,'PaperPositionMode','auto','units','centimeters','Position',[4 4 10 20]);
    set(gcf,'DefaultAxesFontSize',10);
    fontsize = 12;
    linewidth = 1;
    ax1 = axes('position',[.1 .84 .86 .13]);
        plot(td/24/60/60,N_c(xi,:)/1e6,'b','linewidth',linewidth);
        text(0.5,1,'Effective pressure [ MPa ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
        set(gca,'layer','top');
    ax2 = axes('position',[.1 .64 .86 .13]);
        plot(td/24/60/60,Q(xi,:),'b','linewidth',linewidth);
        text(0.5,1,'Discharge [ m${}^3$/s ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
    ax3 = axes('position',[.1 .44 .86 .13]);
        plot(td/24/60/60,S(xi,:),'b','linewidth',linewidth);
        text(0.5,1,'Cross-section [ m${}^2$ ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');
    ax4 = axes('position',[.1 .24 .86 .13]);
        plot(td/24/60/60,M(xi,:),'b','linewidth',linewidth);
        xlabel('Time [ d ]','fontsize',fontsize);
        text(0.5,1,'Melting rate [ m${}^2$/s ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized','fontsize',fontsize,'interpreter','latex');

end

