function f = wyvis_annual_signal(t,t_spr,t_aut,t_per)
% annually periodic function made of two tanh functions centered at tspr and taut with width dt
    ty = 365*24*60*60;
    f = (tanh((mod(t,ty)-t_spr)/t_per)-tanh((mod(t,ty)-t_aut)/t_per))/2;
end