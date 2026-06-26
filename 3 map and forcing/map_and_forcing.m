close all; clc; clear

c_reg{1} = [1, 61, 115] / 255;
c_reg{2} = [58, 174, 140] / 255;
c_reg{3} = [245, 174, 16] / 255;
fs = 9;

gps326 = 'C:\Users\evsalg001\Documents\MATLAB\Zimmer\data\SVP buoys\2025P326_data\2025P326_300534067527860_proc.csv';
gps323 = 'C:\Users\evsalg001\Documents\MATLAB\Zimmer\data\SVP buoys\2025P323_data\2025P323_300534067521880_proc.csv';
gps333 = 'C:\Users\evsalg001\Documents\MATLAB\Zimmer\data\SVP buoys\2025P333_data\2025P333_300534067624300_proc.csv';

sicFolder = "C:\Users\evsalg001\Documents\MATLAB\Zimmer\data\AMSR2_SIC\";
gridFile = sicFolder + "LongitudeLatitudeGrid-n6250-Arctic.hdf";

eraFile = "C:\Users\evsalg001\Documents\MATLAB\Zimmer\data\era5\era5_sr.nc";
surfFile = "C:\Users\evsalg001\Documents\MATLAB\Zimmer\data\Polarstern\PS149_surf_oce.tab";

stationFiles = [
    "C:\Users\evsalg001\Documents\MATLAB\Zimmer\data\radiation stations\2025R29_merged.mat"
    "C:\Users\evsalg001\Documents\MATLAB\Zimmer\data\radiation stations\2025R30_merged.mat"
    "C:\Users\evsalg001\Documents\MATLAB\Zimmer\data\radiation stations\2025R31_merged.mat"
];

stationNames = ["2025R29", "2025R30", "2025R31"];

baseDir = "C:\Users\evsalg001\Documents\MATLAB\Zimmer\data\Polarstern\";
readStation = @(f) readtable(baseDir + f, 'VariableNamingRule','preserve');

T1 = readStation("stations_R1.dat");
T2 = readStation("stations_R2.dat");
T3 = readStation("stations_R3.dat");

lat_1 = T1{:,1}; lon_1 = T1{:,2}; lbl_1 = upper(T1{:,4});
lat_2 = T2{:,1}; lon_2 = T2{:,2}; lbl_2 = upper(T2{:,4});
lat_3 = T3{:,1}; lon_3 = T3{:,2}; lbl_3 = upper(T3{:,4});

Tps = readtable('C:\Users\evsalg001\Documents\MATLAB\Zimmer\data\Polarstern\PS149_mastertrack.tab', 'FileType','text', ...
    'Delimiter','\t', 'HeaderLines',21, 'VariableNamingRule','preserve');

t_ps = datetime(Tps{:,1}, 'InputFormat','yyyy-MM-dd''T''HH:mm');
lat_ps = Tps{:,2};
lon_ps = Tps{:,3};

PS = table(t_ps, lat_ps, lon_ps);
PSd = groupsummary(PS, 't_ps', 'day', {'mean','mean','mean'});

PSd = table( ...
    datetime(string(PSd.day_t_ps)), ...
    PSd.mean_lat_ps, ...
    PSd.mean_lon_ps, ...
    'VariableNames', {'Date','lat','lon'} );

D326 = makeDailyGPS(gps326, datetime(2025,8,21), datetime(2025,7,9),  PSd, "Floe1_P326");
D333 = makeDailyGPS(gps333, datetime(2025,8,24), datetime(2025,7,14), PSd, "Floe2_P333");
D323 = makeDailyGPS(gps323, datetime(2025,8,27), datetime(2025,7,19), PSd, "Floe3_P323");

D = [D326; D333; D323];
D = sortrows(D, {'floe','Date'});
D.Date = dateshift(D.Date, 'start','day');

R = table();

for s = 1:numel(stationFiles)

    S = load(stationFiles(s));
    t_raw = datetime(S.gps.time(:), 'ConvertFrom','datenum');

    Tin = daily_clean_signal( ...
        t_raw, ...
        S.data_f.incom_bb_1, ...
        S.data_f.incom_bb_2, ...
        "MeanFlux");

    Tin.station = repmat(stationNames(s), height(Tin), 1);
    R = [R; Tin];
end

era_time = nc_time_any_units(eraFile, "valid_time");
era_lat = double(ncread(eraFile, "latitude"));
era_lon = double(ncread(eraFile, "longitude"));
era_ssrd = double(ncread(eraFile, "ssrd"));

era_lon = mod(era_lon + 180, 360) - 180;

[era_lon, sortLonIdx] = sort(era_lon);
era_ssrd = era_ssrd(sortLonIdx,:,:);

[era_lat, sortLatIdx] = sort(era_lat);
era_ssrd = era_ssrd(:,sortLatIdx,:);

dt_sec = seconds(median(diff(era_time)));
era_sw = era_ssrd ./ dt_sec;

D.RadStationFlux = NaN(height(D), 1);
D.ERA5Flux = NaN(height(D), 1);
D.MeanFlux = NaN(height(D), 1);
D.RadSource = strings(height(D), 1);
D.RadStation = strings(height(D), 1);

D.RadStation(D.floe == "Floe1_P326") = "2025R29";
D.RadStation(D.floe == "Floe2_P333") = "2025R30";
D.RadStation(D.floe == "Floe3_P323") = "2025R31";

for r = 1:height(D)

    thisDate = dateshift(D.Date(r), 'start','day');
    thisStation = D.RadStation(r);

    idxR = R.Date == thisDate & R.station == thisStation;

    if any(idxR)

        ii = find(idxR, 1);

        D.RadStationFlux(r) = R.MeanFlux(ii);
        D.MeanFlux(r) = R.MeanFlux(ii);
        D.RadSource(r) = thisStation;

    else

        idxEra = dateshift(era_time, 'start','day') == thisDate;

        if any(idxEra)

            era_day = mean(era_sw(:,:,idxEra), 3, 'omitnan');

            lon0 = mod(D.lon(r) + 180, 360) - 180;
            lat0 = D.lat(r);

            D.ERA5Flux(r) = interp2(era_lon, era_lat, era_day', ...
                lon0, lat0, 'linear');

            D.MeanFlux(r) = D.ERA5Flux(r);
            D.RadSource(r) = "ERA5";
        end
    end
end

radius = 2;

lon_sic = double(hdfread(gridFile, 'Longitudes'));
lat_sic = double(hdfread(gridFile, 'Latitudes'));
lon_sic_w = mod(lon_sic + 180, 360) - 180;

F = dir(fullfile(sicFolder, 'asi-AMSR2-n6250-2025*.hdf'));

D.SIC = NaN(height(D), 1);

mapSICsum = NaN(size(lat_sic));
mapSICcount = zeros(size(lat_sic));

for k = 1:length(F)

    fname = F(k).name;
    tok = regexp(fname, 'asi-AMSR2-n6250-(\d{8})-v', 'tokens', 'once');

    if isempty(tok)
        continue
    end

    fileDate = datetime(tok{1}, 'InputFormat','yyyyMMdd');

    sic = double(hdfread(fullfile(F(k).folder, fname), ...
        'ASI Ice Concentration'));
    sic(sic < 0 | sic > 100) = NaN;

    if fileDate >= datetime(2025,7,9) && fileDate <= datetime(2025,8,27)

        valid = isfinite(sic);

        if all(isnan(mapSICsum(:)))
            mapSICsum = zeros(size(sic));
            mapSICcount = zeros(size(sic));
        end

        mapSICsum(valid) = mapSICsum(valid) + sic(valid);
        mapSICcount(valid) = mapSICcount(valid) + 1;
    end

    idx = dateshift(D.Date, 'start','day') == dateshift(fileDate, 'start','day');

    if ~any(idx)
        continue
    end

    rows = find(idx);

    for r = rows'

        lon0 = mod(D.lon(r) + 180, 360) - 180;
        lat0 = D.lat(r);

        dlon = lon_sic_w - lon0;
        dlon = mod(dlon + 180, 360) - 180;

        dist2 = (lat_sic - lat0).^2 + ...
            (dlon .* cosd(lat0)).^2;

        [~, ind] = min(dist2(:));
        [I0, J0] = ind2sub(size(sic), ind);

        Imin = max(1, I0 - radius);
        Imax = min(size(sic,1), I0 + radius);
        Jmin = max(1, J0 - radius);
        Jmax = min(size(sic,2), J0 + radius);

        block = sic(Imin:Imax, Jmin:Jmax);

        D.SIC(r) = mean(block(:), 'omitnan');
    end
end

sic_map = mapSICsum ./ mapSICcount;
sic_map(mapSICcount == 0) = NaN;
sic_map(lat_sic > 87 | lon_sic_w < -35 | lon_sic_w > 50) = NaN;

window_rad = 3;
window_sic = 3;

D.MeanFluxPlot = D.MeanFlux;
D.SIC_smooth = D.SIC;

floeNames = unique(D.floe);

for i = 1:numel(floeNames)

    idx = D.floe == floeNames(i);

    if window_rad > 1
        D.MeanFluxPlot(idx) = movmean(D.MeanFlux(idx), ...
            window_rad, 'omitnan');
    end

    if window_sic > 1
        D.SIC_smooth(idx) = movmean(D.SIC(idx), ...
            window_sic, 'omitnan');
    end
end

idx1 = D.floe == "Floe1_P326";
idx2 = D.floe == "Floe2_P333";
idx3 = D.floe == "Floe3_P323";

eraIdx1 = idx1 & D.RadSource == "ERA5";
eraIdx2 = idx2 & D.RadSource == "ERA5";
eraIdx3 = idx3 & D.RadSource == "ERA5";

lines = readlines(surfFile);
headerLine = find(startsWith(lines,"Date/Time"),1);

opts = detectImportOptions(surfFile, ...
    "FileType","text", ...
    "Delimiter","\t", ...
    "NumHeaderLines",headerLine-1, ...
    "VariableNamingRule","preserve");

Tsurf = readtable(surfFile,opts);

time = datetime(Tsurf.("Date/Time"), ...
    "InputFormat","yyyy-MM-dd'T'HH:mm:ss");

latSurf  = Tsurf.("Latitude");
lonSurf  = Tsurf.("Longitude");
tempSurf = Tsurf.("Temp [°C]");
salSurf  = Tsurf.("Sal");
depthSurf = Tsurf.("Depth water [m]");
tempFlag = Tsurf.("QF water temp");
salFlag  = Tsurf.("QF sal");

p = gsw_p_from_z(-depthSurf, latSurf);
SA = gsw_SA_from_SP(salSurf,p,lonSurf,latSurf);
Tf = gsw_t_freezing(SA,p,0);
deltaT = tempSurf - Tf;

good = ~isnat(time) & ...
       ~isnan(latSurf) & ~isnan(lonSurf) & ...
       ~isnan(tempSurf) & ~isnan(salSurf) & ...
       ~isnan(depthSurf) & ...
       ~isnan(deltaT) & ...
       tempFlag <= 2 & ...
       salFlag <= 2;

time = time(good);
latSurf = latSurf(good);
lonSurf = lonSurf(good);
tempSurf = tempSurf(good);
salSurf = salSurf(good);
Tf = Tf(good);
deltaT = deltaT(good);

visitName = [
    "Visit 1"; "Visit 1"; "Visit 1"
    "Visit 2"; "Visit 2"; "Visit 2"
    "Visit 3"; "Visit 3"; "Visit 3"
    "Visit 4"; "Visit 4"; "Visit 4"
];

floeName = [
    "Floe 1"; "Floe 2"; "Floe 3"
    "Floe 1"; "Floe 2"; "Floe 3"
    "Floe 1 new"; "Floe 2"; "Floe 3"
    "Floe 1 recovery"; "Floe 2"; "Floe 3"
];

startDate = datetime([
    "2025-07-09 06:32"
    "2025-07-14 12:36"
    "2025-07-19 17:45"
    "2025-07-25 16:30"
    "2025-07-30 07:00"
    "2025-08-04 07:28"
    "2025-08-09 13:00"
    "2025-08-12 06:52"
    "2025-08-16 19:55"
    "2025-08-21 16:15"
    "2025-08-23 06:51"
    "2025-08-26 07:10"
], "InputFormat","yyyy-MM-dd HH:mm");

endDate = datetime([
    "2025-07-12 18:25"
    "2025-07-17 17:58"
    "2025-07-22 18:00"
    "2025-07-28 16:40"
    "2025-08-01 20:45"
    "2025-08-06 19:30"
    "2025-08-10 22:07"
    "2025-08-14 17:22"
    "2025-08-18 22:30"
    "2025-08-22 01:00"
    "2025-08-24 18:32"
    "2025-08-27 20:15"
], "InputFormat","yyyy-MM-dd HH:mm");

nStations = numel(visitName);
meanDeltaT = nan(nStations,1);

for i = 1:nStations
    idx = time >= startDate(i) & time <= endDate(i);
    meanDeltaT(i) = mean(deltaT(idx),"omitnan");
end

visitMidDate = startDate + (endDate - startDate)/2;

%%
close all; clc
figure
set(gcf,'Units','inches','Position',[1 3 10.5 5.2])

tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% axMap = nexttile(1,[2 1]);
axMap = nexttile(1);
m_proj('lambert','lons',[-30 40],'lat',[80 85.5]);
hold on

m_pcolor(lon_sic_w, lat_sic, sic_map);
shading flat

load('oslo.mat')
colormap(axMap, oslo)
clim([0 100])

m_gshhs_i('patch',[.7 .7 .7],'edgecolor',[0.3 0.3 0.3],'linewi',0.4);
m_grid('linewi',0.7,'linestyle',':','layer','top','FontSize',fs);

m_line(D.lon(idx1), D.lat(idx1), 'color', c_reg{1}, 'linewi', 2.6);
m_line(D.lon(idx2), D.lat(idx2), 'color', c_reg{2}, 'linewi', 2.6);
m_line(D.lon(idx3), D.lat(idx3), 'color', c_reg{3}, 'linewi', 2.6);

for i = 1:4
    m_line(lon_1(i), lat_1(i), 'marker','o','color','k', ...
        'markerfacecolor',c_reg{1},'linewi',0.8,'markersize',7);

    m_line(lon_2(i), lat_2(i), 'marker','o','color','k', ...
        'markerfacecolor',c_reg{2},'linewi',0.8,'markersize',7);

    m_line(lon_3(i), lat_3(i), 'marker','o','color','k', ...
        'markerfacecolor',c_reg{3},'linewi',0.8,'markersize',7);
end

for i = 1:4
    if i == 3
        p = m_text(lon_1(i)-2.5, lat_1(i)+0.5, lbl_1{i});
        set(p, 'HorizontalAlignment','left', 'FontSize',fs);
    elseif i == 4
        p = m_text(lon_1(i)-1, lat_1(i)-0.2, lbl_1{i});
        set(p, 'HorizontalAlignment','left', 'FontSize',fs);
    else
        p = m_text(lon_1(i), lat_1(i)-0.3, lbl_1{i});
        set(p, 'HorizontalAlignment','center', 'FontSize',fs);
    end
end

for i = 1:4
    p = m_text(lon_2(i)-4.0, lat_2(i)+0.1, lbl_2{i});
    set(p, 'HorizontalAlignment','center', 'FontSize',fs);
end

for i = 1:4
    p = m_text(lon_3(i)+3.5, lat_3(i), lbl_3{i});
    set(p, 'HorizontalAlignment','center', 'FontSize',fs);
end

title('Floe trajectories and sea-ice concentration', ...
    'FontSize',fs, 'FontWeight','normal');

cb = colorbar(axMap);
% cb.Label.String = 'AMSR2 SIC (%)';
cb.Label.String = 'Sea ice concentration (%)';
cb.FontSize = fs;

pos = cb.Position;
pos(1) = pos(1) - 0.04;   % shift left
pos(2) = pos(2) + 0.22;   % shift upward
pos(3) = pos(3) * 0.40;   % make narrower
pos(4) = pos(4) * 0.45;   % make shorter
cb.Position = pos;

axRad = nexttile(2);
hold on

hF1 = plot(D.Date(idx1), D.MeanFluxPlot(idx1), '-', ...
    'Color', c_reg{1}, 'LineWidth', 2);
hF2 = plot(D.Date(idx2), D.MeanFluxPlot(idx2), '-', ...
    'Color', c_reg{2}, 'LineWidth', 2);
hF3 = plot(D.Date(idx3), D.MeanFluxPlot(idx3), '-', ...
    'Color', c_reg{3}, 'LineWidth', 2);

plot(D.Date(eraIdx1), D.MeanFluxPlot(eraIdx1), 'o', ...
    'Color', c_reg{1}, 'MarkerFaceColor', 'w', 'LineWidth', 1.2,'MarkerSize',5)
plot(D.Date(eraIdx2), D.MeanFluxPlot(eraIdx2), 'o', ...
    'Color', c_reg{2}, 'MarkerFaceColor', 'w', 'LineWidth', 1.2,'MarkerSize',5)
plot(D.Date(eraIdx3), D.MeanFluxPlot(eraIdx3), 'o', ...
    'Color', c_reg{3}, 'MarkerFaceColor', 'w', 'LineWidth', 1.2,'MarkerSize',5)

hSmooth = plot(NaN, NaN, '-k', 'LineWidth', 2);
hEra = plot(NaN, NaN, 'ok', 'MarkerFaceColor','w', 'LineWidth',1.2,'MarkerSize',5);

legend([hF1 hF2 hF3 hSmooth hEra], ...
    {'Regime 1','Regime 2','Regime 3', ...
    [num2str(window_rad) '-day moving average'], ...
    'ERA5 gap-fill days'}, ...
    'Location','best', 'FontSize', fs-1,'box','off')

ylabel('Shortwave radiation (W m^{-2})')
ylim([0 250])
box on
axRad.XTickLabel = [];

axSIC = nexttile(4);
hold on

plot(D.Date(idx1), D.SIC_smooth(idx1), '-', ...
    'Color', c_reg{1}, 'LineWidth', 2)
plot(D.Date(idx2), D.SIC_smooth(idx2), '-', ...
    'Color', c_reg{2}, 'LineWidth', 2)
plot(D.Date(idx3), D.SIC_smooth(idx3), '-', ...
    'Color', c_reg{3}, 'LineWidth', 2)

ylabel('Sea ice concentration (%)')
ylim([45 105])
box on

axDT = nexttile(3);
hold on

stationLabel = [
    "1A"; "2A"; "3A"
    "1B"; "2B"; "3B"
    "1C"; "2C"; "3C"
    "1D"; "2D"; "3D"
];

for f = 1:3

    if f == 1
        idx = contains(floeName,'Floe 1');
        floeLabel = 'Floe 1';
    elseif f == 2
        idx = strcmp(floeName,'Floe 2');
        floeLabel = 'Floe 2';
    else
        idx = strcmp(floeName,'Floe 3');
        floeLabel = 'Floe 3';
    end

    x = visitMidDate(idx);
    y = meanDeltaT(idx);
    lbl = stationLabel(idx);

    [x,ord] = sort(x);
    y = y(ord);
     lbl = lbl(ord);

plot(x,y,'-', ...
    'Color',c_reg{f}, ...
    'LineWidth',2)

plot(x,y,'o', ...
    'Color','k', ...
    'MarkerFaceColor',c_reg{f}, ...
    'MarkerSize',7, ...
    'LineWidth',0.8)

    % add labels
for k = 1:numel(x)

    if lbl(k) == "3B"
        dy = -0.04;

    elseif lbl(k) == "2A"
        dy = 0.06;

    elseif lbl(k) == "1D"
        dy = -0.05;

    else
        dy = 0.04;
    end

    text(x(k), y(k)+dy, lbl(k), ...
        'HorizontalAlignment','center', ...
        'FontSize', fs-1)

end
end

ylabel('Surface water T - T_f (°C)')
ylim([0.15 0.85])
box on
datetick(axDT, 'x', 'mmm dd')
xtickangle(axDT, 0)


datetick(axRad, 'x', 'mmm dd')
datetick(axSIC, 'x', 'mmm dd')
xtickangle(axRad, 0)
xtickangle(axSIC, 0)

text(axMap, -0.16, 1.05, '(a)', ...
    'Units','normalized', ...
    'FontSize',10)

text(axRad, -0.12, 1.05, '(b)', ...
    'Units','normalized', ...
    'FontSize',10)

text(axDT, -0.12, 1.05, '(c)', ...
    'Units','normalized', ...
    'FontSize',10)

text(axSIC, -0.12, 1.05, '(d)', ...
    'Units','normalized', ...
    'FontSize',10)

exportgraphics(gcf, 'map_forcing.png', ...
    'Resolution', 300)

%% Helpers

function D = makeDailyGPS(gpsFile, endDate, prependStartDate, PSd, floeName)

G = readtable(gpsFile, 'ReadVariableNames', false);

t = datetime(G.Var1, 'InputFormat','yyyy-MM-dd''T''HH:mm:ss');
lat = G.Var2;
lon = G.Var3;

idx = t <= endDate;
t = t(idx);
lat = lat(idx);
lon = lon(idx);

Tmp = groupsummary(table(t,lat,lon), 't', 'day', ...
    {'mean','mean','mean'}, 'IncludeEmptyGroups',false);

D = table( ...
    forceDatetime(Tmp.day_t), ...
    Tmp.GroupCount, ...
    Tmp.mean_lat, ...
    Tmp.mean_lon, ...
    'VariableNames', {'Date','N','lat','lon'} );

D.floe = repmat(floeName, height(D), 1);

D = prependMissing(D, prependStartDate, PSd, floeName);

end

function D = prependMissing(D, startDate, PSd, floeName)

firstDate = min(D.Date);
missingDates = (dateshift(startDate,'start','day'):days(1):firstDate-days(1))';

if isempty(missingDates)
    return
end

P = table();
P.Date = missingDates;
P.N = NaN(size(missingDates));
P.lat = NaN(size(missingDates));
P.lon = NaN(size(missingDates));
P.floe = repmat(floeName, height(P), 1);

for i = 1:height(P)
    idx = dateshift(PSd.Date, 'start','day') == P.Date(i);
    if any(idx)
        P.lat(i) = PSd.lat(find(idx,1));
        P.lon(i) = PSd.lon(find(idx,1));
    end
end

D = [P; D];

end

function Tdaily = daily_clean_signal(t_raw, rawSignal, cleanSignal, outName)

rawSignal = double(rawSignal(:));
cleanSignal = double(cleanSignal(:));

[~, keepIdx] = match_clean_to_raw(t_raw, rawSignal, cleanSignal);

t = t_raw(keepIdx);
x = cleanSignal(:);

d = dateshift(t, 'start','day');
[Date, ~, idx] = unique(d);

dailyMean = accumarray(idx(:), x(:), [], @mean, NaN);

Tdaily = table(Date, dailyMean, ...
    'VariableNames', {'Date', char(outName)});

end

function [t, keepIdx] = match_clean_to_raw(t_raw, rawSignal, cleanSignal)

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

function t = forceDatetime(x)

if isdatetime(x)
    t = x;
elseif isnumeric(x)
    t = datetime(x, 'ConvertFrom','datenum');
else
    t = datetime(string(x));
end

t = dateshift(t, 'start','day');

end

function t = nc_time_any_units(ncfile, varname)

time_raw = double(ncread(ncfile, varname));
info = ncinfo(ncfile, varname);

units = "";
for i = 1:numel(info.Attributes)
    if strcmpi(info.Attributes(i).Name, 'units')
        units = string(info.Attributes(i).Value);
    end
end

tok = regexp(units, ...
    '(seconds|hours|days) since ([0-9]{4}-[0-9]{2}-[0-9]{2})([ T][0-9]{2}:[0-9]{2}:[0-9]{2})?', ...
    'tokens', 'once');

if isempty(tok)
    t = datetime(time_raw, 'ConvertFrom','posixtime');
else
    unit = tok{1};
    date0 = tok{2};

    if numel(tok) >= 3 && ~isempty(tok{3})
        time0 = strtrim(tok{3});
        t0 = datetime(date0 + " " + time0, ...
            'InputFormat','yyyy-MM-dd HH:mm:ss');
    else
        t0 = datetime(date0, 'InputFormat','yyyy-MM-dd');
    end

    switch unit
        case "seconds"
            t = t0 + seconds(time_raw);
        case "hours"
            t = t0 + hours(time_raw);
        case "days"
            t = t0 + days(time_raw);
    end
end

end