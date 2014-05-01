% Load and clean up Konza data:
clear;
clc;

years = 2006:2012;

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

annual_t = round(100*(1:365*48)/48)/100;

for i = 1:length(years)
    yr = years(i);
    
    % Load the data for that year:
    dtime = ncread(...
        sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',yr),...
        'DTIME') - 1; % -1 because they count have midnight Jan1 = 1
    
    swc1 = ncread(...
        sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',yr),...
        'SWC1'); % in percent
    prec = ncread(...
        sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',yr),...
        'PREC'); % in mm
    ust = ncread(...
        sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',yr),...
        'UST'); % in m/s
    h2o = ncread(...
        sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',yr),...
        'H2O'); % in mmol/mol
    ta = ncread(...
        sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',yr),...
        'TA'); % in C
    h = ncread(...
        sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',yr),...
        'H'); % in W/m^2
    press = ncread(...
        sprintf('Level2/Konza_Prairie/with_gaps/AMF_USKon_%i_L2_WG_V004.nc',yr),...
        'PRESS'); % in kPa
    
    
    
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
    
    
end

SWC1(SWC1 == -9999) = nan;
PREC(PREC == -9999) = nan;
UST(UST == -9999) = nan;
H2O(H2O == -9999) = nan;
TA(TA == -9999) = nan;
H(H == -9999) = nan;
PRESS(PRESS == -9999) = nan;

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

%% Okay, now to work on stability:

stable = zeta > 1;

scatter(SWC1(:),PREC(:));
hold on;
scatter(SWC1(stable),PREC(stable),'.r');



