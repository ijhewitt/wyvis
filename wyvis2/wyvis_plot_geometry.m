function wyvis_plot_geometry(fig,ad)
% plot geometry 
%
% IJH 16 August 2013

    figure(fig); clf;
    set(gcf,'DefaultAxesFontSize',10);
    set(gcf,'Paperpositionmode','auto','units','centimeters','position',[4 4 12 6]);
    linewidth = 2;
    fontsize = 12;
        plot(ad.x/1e3,ad.Z_b,'k','linewidth',linewidth);
        hold on;
        plot(ad.x/1e3,ad.Z_s,'k','linewidth',linewidth);
        if isfield(ad,'x_m')
            xi = wyvis_nearest_gridpoint(ad.x_m,ad.x);
           % plot(ad.x(xi)/1e3*[1 1],[ad.Z_b(xi) ad.Z_s(xi)],'-','color',0.8*[1 1 1],'linewidth',2);
            plot(ad.x(xi)/1e3,ad.Z_b(xi),'bo','MarkerSize',8,'MarkerFaceColor',0.6*[1 1 1]);
        end
        xlabel('Distance [ km ]','fontsize',fontsize,'interpreter','latex');
        ylabel('Elevation [ m ]','fontsize',fontsize,'interpreter','latex');
        set(gca,'layer','top');
        
function xi = wyvis_nearest_gridpoint(x,xd)
% indices of nearest grid points to x
    xi = NaN*x;
    for i = 1:length(x),
        [~,tmp] = min((x(i)-xd).^2);
        xi(i) = tmp;
    end
end
        
end