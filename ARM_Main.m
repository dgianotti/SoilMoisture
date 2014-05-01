% Load and clean up ARM data:
clear;
clc;

years = 2003:2012;

count = reshape(1:(365*length(years)*48),[365*48,length(years)]); % i.e. the time in half hours
t = count/48; % The time in days

YEAR = repmat(years,[365*48,1]);

SWC1 = nan(size(YEAR));
PREC = nan(size(YEAR));
UST = nan(size(YEAR));
H2O = nan(size(YEAR));
TA = nan(size(YEAR));
H = nan(size(YEAR));
PRESS = nan(size(YEAR));
LE = nan(size(YEAR));

annual_t = round(100*(1:365*48)/48)/100;

for i = 1:length(years)
    yr = years(i);
    
    % Load the data for that year:
    dtime = ncread(...
        sprintf('Level2/ARM_SGP_Main/with_gaps/AMF_USARM_%i_L2_WG_V006.nc',yr),...
        'DTIME') - 1; % -1 because they count have midnight Jan1 = 1
    
    swc1 = ncread(...
        sprintf('Level2/ARM_SGP_Main/with_gaps/AMF_USARM_%i_L2_WG_V006.nc',yr),...
        'SWC1'); % in percent
    prec = ncread(...
        sprintf('Level2/ARM_SGP_Main/with_gaps/AMF_USARM_%i_L2_WG_V006.nc',yr),...
        'PREC'); % in mm
    ust = ncread(...
        sprintf('Level2/ARM_SGP_Main/with_gaps/AMF_USARM_%i_L2_WG_V006.nc',yr),...
        'UST'); % in m/s
    h2o = ncread(...
        sprintf('Level2/ARM_SGP_Main/with_gaps/AMF_USARM_%i_L2_WG_V006.nc',yr),...
        'H2O'); % in mmol/mol
    ta = ncread(...
        sprintf('Level2/ARM_SGP_Main/with_gaps/AMF_USARM_%i_L2_WG_V006.nc',yr),...
        'TA'); % in C
    h = ncread(...
        sprintf('Level2/ARM_SGP_Main/with_gaps/AMF_USARM_%i_L2_WG_V006.nc',yr),...
        'H'); % in W/m^2
    press = ncread(...
        sprintf('Level2/ARM_SGP_Main/with_gaps/AMF_USARM_%i_L2_WG_V006.nc',yr),...
        'PRESS'); % in kPa
    le = ncread(...
        sprintf('Level2/ARM_SGP_Main/with_gaps/AMF_USARM_%i_L2_WG_V006.nc',yr),...
        'LE'); % in kPa
    
    
    
    % Round dtime to 2 decimal places:
    dtime = round(dtime*100)/100;
   
    % Put the data into the appropriate matrix:
    SWC1(ismember(annual_t,dtime),i) = swc1(ismember(dtime,annual_t));
    PREC(ismember(annual_t,dtime),i) = prec(ismember(dtime,annual_t));
    UST(ismember(annual_t,dtime),i) = ust(ismember(dtime,annual_t));
    H2O(ismember(annual_t,dtime),i) = h2o(ismember(dtime,annual_t));
    TA(ismember(annual_t,dtime),i) = ta(ismember(dtime,annual_t));
    H(ismember(annual_t,dtime),i) = h(ismember(dtime,annual_t));
    PRESS(ismember(annual_t,dtime),i) = press(ismember(dtime,annual_t));
    LE(ismember(annual_t,dtime),i) = le(ismember(dtime,annual_t));
    
    
end

SWC1(SWC1 == -9999) = nan;
PREC(PREC == -9999) = nan;
UST(UST == -9999) = nan;
H2O(H2O == -9999) = nan;
TA(TA == -9999) = nan;
H(H == -9999) = nan;
PRESS(PRESS == -9999) = nan;
LE(LE == -9999) = nan;

TA = TA + 273.15; % Now iun Kelvins
PRESS = PRESS*1000; % Now in Pa

% Calculate the stability parameter zeta:
z = 3; % air temperature is measured 3m from ground
R = 8.314; % universal gas constant, J/(mol*K)
c_p_air = 1.0035*1000; % J/(kg*K) 
g = 9.8; % m/s^2
partial_press_h2o = H2O/1000;
vk_const = 0.41; % von Karman constant, unitless

air_density = (0.028964*(1-partial_press_h2o) + 0.018016*partial_press_h2o).*PRESS./(R*TA); % in kg/m^3


L_obukhov = (-UST.^3.*air_density*c_p_air.*TA)./(vk_const*g*H);

zeta = z./L_obukhov;
zeta(isinf(zeta)) = nan;

% Remove the daily cycle from LE:
LE_daily_mean = nanmean(reshape(LE,[48,size(LE,2)*365]),2);
LE_daily_anoms = LE - repmat(LE_daily_mean, [365,size(LE,2)]);
LE_daily_annual = LE - repmat( mean(LE,2), [1,size(LE,2)]);

TA_daily_annual = TA - repmat( mean(TA,2), [1,size(TA,2)]);

PREC_daily_annual = PREC - repmat( mean(PREC,2), [1,size(PREC,2)]);

SWC_daily_annual = SWC1 - repmat( mean(SWC1,2), [1,size(SWC1,2)]);

%% Make some plots:

% Daily zeta cycle:
zeta_daily = reshape(zeta,48,size(zeta,2)*365);
x = (24/48):(24/48):24;
quantiles = quantile(zeta_daily,[.05,.5,.95],2);

figure;
plot_between_lines(x',quantiles(:,1),quantiles(:,3),[.75,.75,1])
hold on;
plot(x,quantiles(:,2),'-b','LineWidth',2);
xlim([.5,24]);
ylabel('\zeta');
xlabel('Time of day');
title('Median \zeta with 90% CI');
print(gcf,'-dpng','Zeta_daily_cycle.png');


%%
daytime = abs(t - floor(t) - 0.5) < .25;
%% Some quick correllograms:

%[XC1,lags1] = nancrosscorr(PREC(:),SWC1(:));
%stem(lags1,XC1)

%figure;
%[XC2,lags2] = nancrosscorr(PREC(:),[0;diff(SWC1(:))]);
%stem(lags2,XC2)

figure;
[XC,lags] = nancrosscorr(PREC(:),TA_daily_annual(:));
stem(lags/48,XC)

figure;
[XC,lags] = nancrosscorr(PREC(:),LE_daily_annual(:));
stem(lags/48,XC)

figure;
[XC,lags] = nancrosscorr(PREC(:),SWC1(:));
stem(lags/48,XC)

figure;
[XC,lags] = nancrosscorr(PREC(:),PREC(:));
stem(lags/48,XC)
ylim([0,.1]);

figure;
[XC,lags] = nancrosscorr(SWC_daily_annual(:),LE_daily_annual(:));
stem(lags/48,XC)


%% Okay, now to work on stability:

unstable = zeta > 1;

scatter(SWC1(:),PREC(:));
hold on;
scatter(SWC1(unstable),PREC(unstable),'.r');



