% Time to model some precip!

clear;
clc;

% % % Konza
% Load precip and soil moisture

years = 2006:2012;

SWC1 = [];
PREC= [];
LE = [];
WS = [];
DOY = [];
DTIME = [];

for i = years
     YEAR = ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'YEAR');
%     ZL = ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'ZL');
% 
%     % Water vapor content
%     H2O = ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'H2O');
%     H2O(H2O == -9999) = nan;
%     
%     Rgl = ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'Rgl');
%     RglOut = ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'RglOut');
%     RgOut = ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'RgOut');
%     
%     % Net radiation
%     Rn = ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'Rn');
%     Rn(Rn==-9999) = nan;
    
    SWC1 = [SWC1; ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'SWC1')];
    
%     SWC2 = ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'SWC2');
%     SWC2(SWC2 == -9999) = nan;
%     
%     PRESS = ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'PRESS');
%     RH = ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'RH');
    
    % Precipitation
    PREC = [PREC; ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'PREC')];
    
%     TS1 = ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'TS1');
%     TSdepth1 = ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'TSdepth1');
%     TS2 = ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'TS2');
%     TSdepth2 = ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'TSdepth2');
%     FG = ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'FG');
    
    % Latent heat flux
    LE = [LE; ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'LE')];
    
    WS = [WS; ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'WS')];
    
%     WD = ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'WD');
%     TA = ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'TA');
    
    % u*
%     UST = ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'UST');
%     UST(UST == -9999) = nan;
    
%     HRMIN = ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'HRMIN');
DOY = [DOY; ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'DOY')];
    DTIME = [DTIME; ncread(sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',i),'DTIME')];
    

end

SWC1(SWC1 == -9999) = nan;
WS(WS == -9999) = nan;
LE(LE == -9999) = nan;
PREC(PREC == -9999) = nan;
    

data = PREC;
scaled = (data - nanmean(data)) / nanstd(data);
scaled(isnan(data)) = 0;
[corr,lags] = xcorr(scaled,'coeff');
stem(lags,corr);

data = SWC1;
scaled = (data - nanmean(data)) / nanstd(data);
scaled(isnan(data)) = 0;
[corr,lags] = xcorr(scaled,'coeff');
stem(lags,corr);



[XC1,lags1] = nancrosscorr(PREC,SWC1);
stem(lags1,XC1)

figure;
[XC2,lags2] = nancrosscorr(PREC,[0;diff(SWC1)]);
stem(lags2,XC2)

[a,b] = max(XC2)
lags2(b)


%% Need to make cleaner data set:

