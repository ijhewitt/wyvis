clear;

fn = 'KAN_M_day_v03.txt';
lat = 67.0670; lon = -48.8355; el = 1270;
M = dlmread(fn,'',1,0);
% Y = M(:,1); % year
% DOY = M(:,4); % day of year
D = M(:,5); % day of century
T = M(:,7); T(T<-100) = NaN; % temp (C)

% 1 of Jan
DtoY = [
2923    2008;
3289    2009;
3654    2010;
4019   2011;
4384    2012;
4750    2013;
5115    2014;
5480    2015;
5845    2016;
6211    2017;
6576    2018;
6941    2019];

figure(1); clf;
    plot(D,0*D.^0,'k--'); hold on; 
    plot(D,T,'.-');
    set(gca,'XTick',DtoY(:,1),'XTickLabel',DtoY(:,2));
    xlabel('Time');
    ylabel('T [C]');

N = 7; Tsmooth = filter(1/N*ones(1,N),1,T); % moving average

% figure(1);
%     plot(D,Tsmooth); 
%     xlim([5480 6211]);
    

