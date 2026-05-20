function wyvis_movie_u(fig,ad,ud,td)
% wyvis_movie_u(fig,ad,ud,td)
% plot movie of u from solution structure ud
    td2 = 24*60*60;
    for i = 1:length(td), 
        wyvis_plot_u(fig,ad,ud(i)); 
        if i==1, axis0 = axis; else axis(axis0); end
        text(0.05,0.1,['t = ',num2str(round(td(i)/td2)),'d'],'units','normalized'); 
        drawnow; shg;
    end
end