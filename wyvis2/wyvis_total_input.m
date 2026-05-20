function in = wyvis_total_input(td,ad)
% integrates all inputs with respect to x 
    in = cumtrapz(ad.x,ad.M_in(td),1) ...
        + cumtrapz(ad.x.*ad.W,ad.m_in(td),1);
    if length(ad.Q_m(0))>1,
        in = in + cumsum(ad.Q_in(td),1); 
    end
    in = in(end,:);
end