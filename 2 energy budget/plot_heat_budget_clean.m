%% plot_heat_budget_clean.m
%
% Clean script that only produces the final heat-budget figure from the original script:
%   (a) basal melt heat flux / OHF
%   (b) open-water absorbed shortwave flux
%   (c) transmitted shortwave flux through sea ice
%
% Required input files:
%   - Floe GPS CSV files
%   - AMSR2 SIC HDF files and LongitudeLatitudeGrid-n6250-Arctic.hdf
%   - ERA5 SSRD NetCDF file
%   - Radiation station merged MAT files
%   - SIMBA NetCDF files
%   - PS149_mastertrack.tab
%
% Update paths in the USER SETTINGS block if needed.

clear; clc; close all

%% USER SETTINGS

% --- Floe GPS files
fileGPS326 = 'C:\Users\evsalg001\Documents\MATLAB\Zimmer\data\SVP buoys\2025P326_data\2025P326_300534067527860_proc.csv';
fileGPS333 = 'C:\Users\evsalg001\Documents\MATLAB\Zimmer\data\SVP buoys\2025P333_data\2025P333_300534067624300_proc.csv';
fileGPS323 = 'C:\Users\evsalg001\Documents\MATLAB\Zimmer\data\SVP buoys\2025P323_data\2025P323_300534067521880_proc.csv';

% --- Polarstern master track
filePS = 'C:\Users\evsalg001\Documents\MATLAB\Zimmer\data\Polarstern\PS149_mastertrack.tab';

% --- AMSR2 SIC
sicFolder = 'C:\Users\evsalg001\Documents\MATLAB\Zimmer\data\AMSR2_SIC';
gridFile  = fullfile(sicFolder, 'LongitudeLatitudeGrid-n6250-Arctic.hdf');
sicRadiusPixels = 2;

% --- ERA5 shortwave radiation
fileERA5 = 'C:\Users\evsalg001\Documents\MATLAB\Zimmer\data\era5\era5_sr.nc';

% --- Radiation stations
stationFiles = [
    "C:\Users\evsalg001\Documents\MATLAB\Zimmer\data\radiation stations\2025R29_merged.mat"
    "C:\Users\evsalg001\Documents\MATLAB\Zimmer\data\radiation stations\2025R30_merged.mat"
    "C:\Users\evsalg001\Documents\MATLAB\Zimmer\data\radiation stations\2025R31_merged.mat"
];

stationNames = ["2025R29"; "2025R30"; "2025R31"];

% --- SIMBA data
simbaDir = 'C:\Users\evsalg001\Documents\MATLAB\Zimmer\data\SIMBA';
buoyIDs = {'T143','T144','T145','T135','T136'};

% --- Physics and plotting
alpha_ow = 0.066;        % open-water albedo
Lvol = 260e6;            % effective volumetric latent heat [J m^-3]
smoothWindow = 7;        % smoothing window for final plot [days]
exportDPI = 300;
outputFile = 'heat_budget_components.png';

% Regime colors
c_reg{1} = [1, 61, 115] / 255;
c_reg{2} = [58, 174, 140] / 255;
c_reg{3} = [245, 174, 16] / 255;

% Floe labels and dates
floeNames = ["Floe1_P326"; "Floe2_P333"; "Floe3_P323"];

% GPS end dates used in the original script
endDate326 = datetime(2025,8,21);
endDate333 = datetime(2025,8,24);
endDate323 = datetime(2025,8,27);

% Dates before floe GPS are filled with daily Polarstern positions
startDate326 = datetime(2025,7,9);
startDate333 = datetime(2025,7,14);
startDate323 = datetime(2025,7,19);

%% READ POLARSTERN DAILY POSITIONS

Tps = readtable(filePS, 'FileType','text', ...
    'Delimiter','\t', 'HeaderLines',21, 'VariableNamingRule','preserve');

t_ps = datetime(Tps{:,1}, 'InputFormat','yyyy-MM-dd''T''HH:mm');
lat_ps = Tps{:,2};
lon_ps = Tps{:,3};

PS = table(t_ps, lat_ps, lon_ps, ...
    'VariableNames', {'time','lat','lon'});

PSd0 = groupsummary(PS, 'time', 'day', {'mean','mean'});
PSd = table( ...
    forceDatetime(PSd0.day_time), ...
    PSd0.mean_lat, ...
    PSd0.mean_lon, ...
    'VariableNames', {'Date','lat','lon'});

%% READ FLOE GPS AND MAKE DAILY TABLE

D326 = readFloeGPS(fileGPS326, endDate326, startDate326, PSd, "Floe1_P326");
D333 = readFloeGPS(fileGPS333, endDate333, startDate333, PSd, "Floe2_P333");
D323 = readFloeGPS(fileGPS323, endDate323, startDate323, PSd, "Floe3_P323");

D = [D326; D333; D323];
D = sortrows(D, {'floe','Date'});

%% DAILY AMSR2 SIC AT FLOE POSITIONS

lon_sic = double(hdfread(gridFile, 'Longitudes'));
lat_sic = double(hdfread(gridFile, 'Latitudes'));
lon_sic = wrapTo180local(lon_sic);

F = dir(fullfile(sicFolder, 'asi-AMSR2-n6250-2025*.hdf'));
D.SIC = NaN(height(D), 1);

for k = 1:numel(F)

    fname = F(k).name;
    tok = regexp(fname, 'asi-AMSR2-n6250-(\d{8})-v', 'tokens', 'once');

    if isempty(tok)
        continue
    end

    fileDate = datetime(tok{1}, 'InputFormat','yyyyMMdd');
    idxDate = dateshift(D.Date, 'start','day') == dateshift(fileDate, 'start','day');

    if ~any(idxDate)
        continue
    end

    sic = double(hdfread(fullfile(F(k).folder, fname), 'ASI Ice Concentration'));
    sic(sic < 0 | sic > 100) = NaN;

    rows = find(idxDate);

    for r = rows'

        lon0 = wrapTo180local(D.lon(r));
        lat0 = D.lat(r);

        dlon = wrapTo180local(lon_sic - lon0);
        dist2 = (lat_sic - lat0).^2 + (dlon .* cosd(lat0)).^2;

        [~, ind] = min(dist2(:));
        [I0, J0] = ind2sub(size(sic), ind);

        Imin = max(1, I0 - sicRadiusPixels);
        Imax = min(size(sic,1), I0 + sicRadiusPixels);
        Jmin = max(1, J0 - sicRadiusPixels);
        Jmax = min(size(sic,2), J0 + sicRadiusPixels);

        block = sic(Imin:Imax, Jmin:Jmax);
        D.SIC(r) = mean(block(:), 'omitnan');
    end
end

D.SIC_smooth = NaN(height(D), 1);
for i = 1:numel(floeNames)
    idx = D.floe == floeNames(i);
    D.SIC_smooth(idx) = movmean(D.SIC(idx), 5, 'omitnan');
end

%% DAILY INCOMING SHORTWAVE FROM RADIATION STATIONS

R = table();

for s = 1:numel(stationFiles)

    S = load(stationFiles(s));

    t_raw = datetime(S.gps.time(:), 'ConvertFrom','datenum');
    rad_raw = double(S.data_f.incom_bb_1(:));
    rad_clean = double(S.data_f.incom_bb_2(:));

    [t, keepIdx] = matchCleanToRaw(t_raw, rad_raw, rad_clean);
    rad = rad_clean(:);

    lat = double(S.gps.LAT(:));
    lon = wrapTo180local(double(S.gps.LON(:)));

    if isscalar(lat)
        lat = repmat(lat, size(t));
    else
        lat = lat(keepIdx);
    end

    if isscalar(lon)
        lon = repmat(lon, size(t));
    else
        lon = lon(keepIdx);
    end

    Tst = table(t, lat, lon, rad, ...
        'VariableNames', {'time','lat','lon','rad'});

    Td = groupsummary(Tst, 'time', 'day', {'mean','mean','mean'});

    Rtmp = table();
    Rtmp.Date = forceDatetime(Td.day_time);
    Rtmp.N = Td.GroupCount;
    Rtmp.lat = Td.mean_lat;
    Rtmp.lon = Td.mean_lon;
    Rtmp.MeanFlux = Td.mean_rad;
    Rtmp.station = repmat(stationNames(s), height(Rtmp), 1);

    R = [R; Rtmp]; %#ok<AGROW>
end

%% ERA5 DAILY SHORTWAVE FOR GAP FILLING

era_time = ncTime(fileERA5, 'valid_time');
era_lat = double(ncread(fileERA5, 'latitude'));
era_lon = wrapTo180local(double(ncread(fileERA5, 'longitude')));
era_ssrd = double(ncread(fileERA5, 'ssrd'));

[era_lon, sortLonIdx] = sort(era_lon);
era_ssrd = era_ssrd(sortLonIdx,:,:);

[era_lat, sortLatIdx] = sort(era_lat);
era_ssrd = era_ssrd(:,sortLatIdx,:);

dt_sec = seconds(median(diff(era_time)));
era_sw = era_ssrd ./ dt_sec;

%% MERGE RADIATION, SIC, OPEN-WATER SW, AND TRANSMITTED SW

D2 = D;
D2.Date = dateshift(D2.Date, 'start','day');
D2.MeanFlux = NaN(height(D2),1);
D2.RadSource = strings(height(D2),1);
D2.RadStation = strings(height(D2),1);

D2.RadStation(D2.floe == "Floe1_P326") = "2025R29";
D2.RadStation(D2.floe == "Floe2_P333") = "2025R30";
D2.RadStation(D2.floe == "Floe3_P323") = "2025R31";

for r = 1:height(D2)

    thisDate = dateshift(D2.Date(r), 'start','day');
    thisStation = D2.RadStation(r);

    idxR = R.Date == thisDate & R.station == thisStation;

    if any(idxR)
        ii = find(idxR, 1);
        D2.MeanFlux(r) = R.MeanFlux(ii);
        D2.RadSource(r) = thisStation;
    else
        idxEra = dateshift(era_time, 'start','day') == thisDate;

        if any(idxEra)
            era_day = mean(era_sw(:,:,idxEra), 3, 'omitnan');

            lon0 = wrapTo180local(D2.lon(r));
            lat0 = D2.lat(r);

            D2.MeanFlux(r) = interp2( ...
                era_lon, era_lat, era_day', ...
                lon0, lat0, 'linear');

            D2.RadSource(r) = "ERA5";
        end
    end
end

% Open-water absorbed shortwave radiation
D2.f_ow = 1 - D2.SIC_smooth ./ 100;
D2.F_into_ow = D2.MeanFlux .* D2.f_ow .* (1 - alpha_ow);

% Regime labels
D2.regime = NaN(height(D2),1);
D2.regime(D2.floe == "Floe1_P326") = 1;
D2.regime(D2.floe == "Floe2_P333") = 2;
D2.regime(D2.floe == "Floe3_P323") = 3;

% ROV-derived broadband transmittance at visit dates.
% Regime 3 values use the corrected transmittance series.
Trov = table();
Trov.Date = [
    datetime(2025,7,11,12,0,0)
    datetime(2025,7,27,12,0,0)
    datetime(2025,8,9,12,0,0)
    datetime(2025,7,16,12,0,0)
    datetime(2025,7,31,12,0,0)
    datetime(2025,8,12,12,0,0)
    datetime(2025,8,23,12,0,0)
    datetime(2025,7,22,12,0,0)
    datetime(2025,8,4,12,0,0)
    datetime(2025,8,17,12,0,0)
    datetime(2025,8,26,12,0,0)
];

Trov.regime = [
    1; 1; 1
    2; 2; 2; 2
    3; 3; 3; 3
];

Trov.Transmission = [
    0.2082
    0.1018
    0.0891
    0.0662
    0.0986
    0.1169
    0.1019
    0.0538
    0.0601
    0.0552
    0.0562
];

D2.Transmission = NaN(height(D2),1);

for reg = 1:3

    idxD = D2.regime == reg;
    idxT = Trov.regime == reg;

    D2.Transmission(idxD) = interp1( ...
        datenum(Trov.Date(idxT)), ...
        Trov.Transmission(idxT), ...
        datenum(D2.Date(idxD)), ...
        'linear', 'extrap');
end

D2.Transmission(D2.Transmission < 0) = NaN;

sic_frac = D2.SIC_smooth ./ 100;
D2.F_transm = D2.MeanFlux .* sic_frac .* D2.Transmission;

%% SIMBA BASAL MELT HEAT FLUX

regimeMap = containers.Map( ...
    {'T143','T144','T145','T135','T136'}, ...
    [1      2      3      2      3]);

B = table();

for b = 1:numel(buoyIDs)

    buoyID = buoyIDs{b};
    ncFile = fullfile(simbaDir, sprintf('2025%s.nc', buoyID));

    t = epochDaysToDatetime(ncread(ncFile, 'time_temperature'));
    iceBottom = double(ncread(ncFile, 'ice_water_interface_temperature'));

    t = t(:);
    iceBottom = iceBottom(:);

    dt_days = days(diff(t));
    tMid = t(1:end-1) + diff(t)/2;

    basalMelt = -diff(iceBottom) ./ dt_days;    % m day^-1
    basalFlux = Lvol .* basalMelt ./ 86400;     % W m^-2

    good = isfinite(basalFlux);

    Tm = table( ...
        dateshift(tMid(good), 'start','day'), ...
        basalMelt(good), ...
        basalFlux(good), ...
        'VariableNames', {'Date','BasalMelt','BasalFlux'});

    Td = groupsummary(Tm, 'Date', 'mean', {'BasalMelt','BasalFlux'});

    tmp = table();
    tmp.Date = Td.Date;
    tmp.BasalMelt = Td.mean_BasalMelt;
    tmp.BasalFlux = Td.mean_BasalFlux;
    tmp.buoy = repmat(string(buoyID), height(tmp), 1);
    tmp.regime = repmat(regimeMap(buoyID), height(tmp), 1);

    B = [B; tmp]; %#ok<AGROW>
end

% Keep only dates where radiative fluxes exist
fluxDates = unique(D2.Date(~isnan(D2.F_into_ow) | ~isnan(D2.F_transm)));
B = B(ismember(B.Date, fluxDates), :);

% Mean by regime and date
Breg = groupsummary(B, {'Date','regime'}, 'mean', {'BasalMelt','BasalFlux'});

%% DAILY HEAT BUDGET BY REGIME

Budget = table();

for reg = 1:3

    thisFloe = floeNames(reg);

    idxF = D2.floe == thisFloe;
    idxB = Breg.regime == reg;

    Tsw = table();
    Tsw.Date = D2.Date(idxF);
    Tsw.F_ow = D2.F_into_ow(idxF);
    Tsw.F_transm = D2.F_transm(idxF);
    Tsw.F_sw_total = Tsw.F_ow + Tsw.F_transm;

    Tohf = table();
    Tohf.Date = Breg.Date(idxB);
    Tohf.F_ohf = Breg.mean_BasalFlux(idxB);

    Tb = innerjoin(Tohf, Tsw, 'Keys','Date');
    Tb.regime = repmat(reg, height(Tb), 1);

    Budget = [Budget; Tb]; %#ok<AGROW>
end

% Smoothed terms for final figure
Budget.F_ohf_smooth = NaN(height(Budget),1);
Budget.F_ow_smooth = NaN(height(Budget),1);
Budget.F_transm_smooth = NaN(height(Budget),1);

for reg = 1:3
    idx = Budget.regime == reg;

    Budget.F_ohf_smooth(idx) = movmean(Budget.F_ohf(idx), ...
        smoothWindow, 'omitnan');

    Budget.F_ow_smooth(idx) = movmean(Budget.F_ow(idx), ...
        smoothWindow, 'omitnan');

    Budget.F_transm_smooth(idx) = movmean(Budget.F_transm(idx), ...
        smoothWindow, 'omitnan');
end

%% FINAL FIGURE: 3-PANEL HEAT BUDGET

figure
set(gcf, 'Units','inches', 'Position',[3 3 8 6])

tl = tiledlayout(3,1, 'TileSpacing','compact', 'Padding','compact');

% Panel a: basal OHF
ax1 = nexttile;
hold on; box on
for reg = 1:3
    idx = Budget.regime == reg;
    plot(Budget.Date(idx), Budget.F_ohf_smooth(idx), '-', ...
        'Color', c_reg{reg}, 'LineWidth', 2.5)
end
ylabel({'Basal melt heat flux'; '(W m^{-2})'}, 'FontSize', 9)
ylim([0 200])
legend({'Regime 1','Regime 2','Regime 3'}, ...
    'Location','north', 'Box','off', 'FontSize',8, 'Orientation','horizontal')
addMeanText(ax1, Budget, 'F_ohf', 'Mean OHF', c_reg, 'right')

% Panel b: open-water SW
ax2 = nexttile;
hold on; box on
for reg = 1:3
    idx = Budget.regime == reg;
    plot(Budget.Date(idx), Budget.F_ow_smooth(idx), '-', ...
        'Color', c_reg{reg}, 'LineWidth', 2.5)
end
ylabel({'Open-water SW'; '(W m^{-2})'}, 'FontSize', 9)
ylim([0 40])
addMeanText(ax2, Budget, 'F_ow', 'Mean OW SW', c_reg, 'left')

% Panel c: transmitted SW
ax3 = nexttile;
hold on; box on
for reg = 1:3
    idx = Budget.regime == reg;
    plot(Budget.Date(idx), Budget.F_transm_smooth(idx), '-', ...
        'Color', c_reg{reg}, 'LineWidth', 2.5)
end
ylabel({'Transmitted SW'; '(W m^{-2})'}, 'FontSize', 9)
ylim([0 40])
addMeanText(ax3, Budget, 'F_transm', 'Mean transmitted SW', c_reg, 'right')

% Axes formatting

datetick(ax1,'x','mmm dd')
datetick(ax2,'x','mmm dd')
datetick(ax3,'x','mmm dd')

ax1.XTickLabel = [];
ax2.XTickLabel = [];

linkaxes([ax1 ax2 ax3], 'x')

text(ax1,-0.07,1.03,'(a)', ...
    'Units','normalized', ...
    'FontSize',9);

text(ax2,-0.07,1.03,'(b)', ...
    'Units','normalized', ...
    'FontSize',9);

text(ax3,-0.07,1.03,'(c)', ...
    'Units','normalized', ...
    'FontSize',9);

exportgraphics(gcf, outputFile, 'Resolution', exportDPI)

%% PRINT MEAN VALUES TO COMMAND WINDOW

fprintf('\nMean fluxes over matched heat-budget period:\n')
fprintf('Regime   OHF       OW SW     Transm SW\n')
fprintf('         W m^-2    W m^-2    W m^-2\n')
for reg = 1:3
    idx = Budget.regime == reg;
    fprintf('R%d       %6.1f    %6.1f    %6.1f\n', ...
        reg, ...
        mean(Budget.F_ohf(idx), 'omitnan'), ...
        mean(Budget.F_ow(idx), 'omitnan'), ...
        mean(Budget.F_transm(idx), 'omitnan'))
end
fprintf('\nSaved figure: %s\n', outputFile)

%% LOCAL FUNCTIONS

function Dnew = readFloeGPS(fileGPS, endDate, startDate, PSd, floeName)

G = readtable(fileGPS, 'ReadVariableNames', false);

t = datetime(G.Var1, 'InputFormat','yyyy-MM-dd''T''HH:mm:ss');
lat = G.Var2;
lon = G.Var3;

idx = t <= endDate;
t = t(idx);
lat = lat(idx);
lon = lon(idx);

Tin = table(t, lat, lon, 'VariableNames', {'t','lat','lon'});
Td = groupsummary(Tin, 't', 'day', {'mean','mean'});

Dfloe = table( ...
    forceDatetime(Td.day_t), ...
    Td.GroupCount, ...
    Td.mean_lat, ...
    Td.mean_lon, ...
    'VariableNames', {'Date','N','lat','lon'});

Dnew = prependMissing(Dfloe, startDate, PSd, floeName);

end

function Dnew = prependMissing(Dfloe, startDate, PSd, floeName)

allDates = (startDate : days(1) : max(Dfloe.Date))';

Dnew = table( ...
    allDates, ...
    NaN(size(allDates)), ...
    NaN(size(allDates)), ...
    NaN(size(allDates)), ...
    'VariableNames', {'Date','N','lat','lon'});

[~, ia, ib] = intersect(Dnew.Date, Dfloe.Date);
Dnew.N(ia) = Dfloe.N(ib);
Dnew.lat(ia) = Dfloe.lat(ib);
Dnew.lon(ia) = Dfloe.lon(ib);

missing = isnan(Dnew.lat);
[~, ia2, ib2] = intersect(Dnew.Date(missing), PSd.Date);

idxMissing = find(missing);
Dnew.lat(idxMissing(ia2)) = PSd.lat(ib2);
Dnew.lon(idxMissing(ia2)) = PSd.lon(ib2);
Dnew.N(idxMissing(ia2)) = 1;

Dnew.floe = repmat(floeName, height(Dnew), 1);

end

function dt = forceDatetime(x)

if isdatetime(x)
    dt = x;
elseif iscell(x)
    dt = vertcat(x{:});
elseif iscategorical(x)
    dt = datetime(string(x));
elseif isstring(x) || ischar(x)
    dt = datetime(x);
else
    error('Unsupported date format in groupsummary output')
end

dt = dateshift(dt, 'start','day');

end

function t = ncTime(ncfile, varname)

x = double(ncread(ncfile, varname));
u = string(ncreadatt(ncfile, varname, 'units'));

ref = erase(extractAfter(u, 'since '), '.0');

if startsWith(u, 'hours')
    t = datetime(ref, 'InputFormat','yyyy-MM-dd HH:mm:ss') + hours(x);
elseif startsWith(u, 'seconds')
    t = datetime(ref, 'InputFormat','yyyy-MM-dd HH:mm:ss') + seconds(x);
elseif startsWith(u, 'days')
    t = datetime(ref, 'InputFormat','yyyy-MM-dd HH:mm:ss') + days(x);
else
    error('Unsupported NetCDF time units: %s', u)
end

t = t(:);

end

function t = epochDaysToDatetime(timeDays)

timeDays = double(timeDays(:));
t = datetime(1970,1,1,0,0,0) + days(timeDays);

end

function [t, keepIdx] = matchCleanToRaw(t_raw, rawSignal, cleanSignal)

keepIdx = NaN(size(cleanSignal));
j0 = 1;

for j = 1:numel(cleanSignal)

    ii = find(abs(rawSignal(j0:end) - cleanSignal(j)) < 1e-10, 1, 'first');

    if isempty(ii)
        error('Could not match clean signal value %d to raw signal.', j)
    end

    keepIdx(j) = j0 + ii - 1;
    j0 = keepIdx(j) + 1;
end

t = t_raw(keepIdx);

end

function lon = wrapTo180local(lon)

lon = mod(lon + 180, 360) - 180;

end

function addMeanText(ax, Budget, varName, titleText, c_reg, side)

m1 = mean(Budget.(varName)(Budget.regime == 1), 'omitnan');
m2 = mean(Budget.(varName)(Budget.regime == 2), 'omitnan');
m3 = mean(Budget.(varName)(Budget.regime == 3), 'omitnan');

txt = sprintf('%s\nR1: %.1f W m^{-2}\nR2: %.1f W m^{-2}\nR3: %.1f W m^{-2}', ...
    titleText, m1, m2, m3);

if strcmpi(side, 'right')
    x = 0.98;
    hAlign = 'right';
else
    x = 0.02;
    hAlign = 'left';
end

text(ax, x, 0.95, txt, ...
    'Units','normalized', ...
    'HorizontalAlignment', hAlign, ...
    'VerticalAlignment','top', ...
    'FontSize', 8)

% Add small colored regime markers near the mean text for visual consistency.
% Kept deliberately simple to avoid clutter.
unusedColors = c_reg; %#ok<NASGU>

end
