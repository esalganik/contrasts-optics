% plot_orthomosaic.m
%
% Visualization of orthomosaics.
%
% Developed for:
% "Contrasting Optical Properties of Different Arctic Sea Ice Regimes"

clear; clc; close all

% Path to the orthomosaic GeoTIFF.
% Change this path if the file is stored in a different location.

filename = 'C:\Users\evsalg001\Downloads\1b_georeferenced.tif';

% Image subset to plot

xlim_patch = [-100 250];     % x-range of displayed orthomosaic [m]
ylim_patch = [0 250];        % y-range of displayed orthomosaic [m]

% Figure export

fig_width_in  = 6;           % exported figure width [in]
fig_height_in = 4;           % exported figure height [in]

export_dpi = 300;            % export resolution [dpi]

output_file = 'orthomosaic_1B.tif';

% Albedo transect

y_profile = 170;             % y-position of albedo transect [m]

x_dark = 95;                 % position of open-water reference point [m]
albedo_dark = 0.07;          % albedo assigned to open water

x_bright = 145;              % position of bright ice reference point [m]
albedo_bright = 0.70;        % albedo assigned to bright ice

% Albedo marker display

marker_spacing_m = 3;        % spacing between plotted albedo markers [m]

%% Import image patch

ginfo = georasterinfo(filename);
R = ginfo.RasterReference;

col1 = floor((xlim_patch(1) - R.XWorldLimits(1)) / R.CellExtentInWorldX) + 1;
col2 = ceil((xlim_patch(2) - R.XWorldLimits(1)) / R.CellExtentInWorldX);

row1 = floor((R.YWorldLimits(2) - ylim_patch(2)) / R.CellExtentInWorldY) + 1;
row2 = ceil((R.YWorldLimits(2) - ylim_patch(1)) / R.CellExtentInWorldY);

col1 = max(1,col1);
col2 = min(R.RasterSize(2),col2);

row1 = max(1,row1);
row2 = min(R.RasterSize(1),row2);

target_width_px  = fig_width_in  * export_dpi;
target_height_px = fig_height_in * export_dpi;

ncols_full = col2 - col1 + 1;
nrows_full = row2 - row1 + 1;

decim_x = ceil(ncols_full / target_width_px);
decim_y = ceil(nrows_full / target_height_px);
decim = max([1 decim_x decim_y]);

cols = col1:decim:col2;
rows = row1:decim:row2;

fprintf('Original patch: %d x %d pixels\n',ncols_full,nrows_full)
fprintf('Imported patch: %d x %d pixels\n',length(cols),length(rows))
fprintf('Decimation factor: %d\n',decim)

A = imread(filename,...
    'PixelRegion',{[row1 decim row2],[col1 decim col2]});

RGB = A(:,:,1:3);
mask = A(:,:,4) > 0;

x = R.XWorldLimits(1) + ...
    (cols - 0.5) * R.CellExtentInWorldX;

y = R.YWorldLimits(2) - ...
    (rows - 0.5) * R.CellExtentInWorldY;

%% Calculate albedo transect

[~,row_profile] = min(abs(y - y_profile));

brightness = squeeze(mean(double(RGB(row_profile,:,:)),3));
brightness(~mask(row_profile,:)) = NaN;

[~,i_dark] = min(abs(x - x_dark));
[~,i_bright] = min(abs(x - x_bright));

B_dark = brightness(i_dark);
B_bright = brightness(i_bright);

if B_bright == B_dark
    error('Calibration points have identical brightness.')
end

albedo = albedo_dark + ...
    (brightness - B_dark) .* ...
    (albedo_bright - albedo_dark) ./ ...
    (B_bright - B_dark);

albedo(albedo < 0) = 0;
albedo(albedo > 1) = 1;

%% Plot

fig = figure('Color','w');
fig.Units = 'inches';
fig.Position = [1 4 fig_width_in fig_height_in];

h = image(x,y,RGB);
set(h,'AlphaData',mask)

axis image
set(gca,'YDir','normal')
hold on

marker_step = max(1,...
    round(marker_spacing_m/(decim*R.CellExtentInWorldX)));

idx = 1:marker_step:length(x);

scatter(x(idx), ...
        y_profile*ones(size(idx)), ...
        12, ...
        albedo(idx), ...
        'filled', ...
        'MarkerEdgeColor','none')

% Optional: if buda.mat (Crameri colormap) is available in the MATLAB path,
% it will be used. Otherwise MATLAB's default parula colormap is used.
if exist('buda.mat','file')
    S = load('buda.mat');
    colormap(S.buda)
else
    colormap(parula)
end

clim([0 1])

cb = colorbar;
cb.Label.String = 'Albedo';

plot(146.2,148.6,'r+','MarkerSize',4,'LineWidth',1)
text(146.2,148.6,'  ROV hole','Color','r','FontSize',7)

plot(176.6,132.3,'r+','MarkerSize',4,'LineWidth',1)
text(176.6,132.3,'  Radiation Station','Color','r','FontSize',7)

xlim(xlim_patch)
ylim(ylim_patch)

xlabel('x (m)')
ylabel('y (m)')
title('Floe 1B orthomosaic with albedo transect',...
      'FontWeight','normal')

% Export

exportgraphics(fig,output_file,'Resolution',export_dpi)