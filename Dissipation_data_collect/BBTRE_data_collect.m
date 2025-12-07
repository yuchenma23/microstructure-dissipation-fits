clear;
clc;
%% --- USER OPTIONS ------------------------------------------------------
base_path_96 = '../../../../Data/BBTRE/BBTRE_data/BBTRE96/bbtre96_';
base_path_97 = '../../../../Data/BBTRE/BBTRE_data/BBTRE97/bbtre97_';

lon_min = -27; lon_max = -10;
lat_min = -30; lat_max =  -6;

edges = [0 1000 2000];      % HAB bins: 0–1000 m, 1000–2000 m
useRegionFilter = true;

% two output vectors of dissipation
epsilons_bbtre_0_1000   = [];   % HAB in [0, 1000) m
epsilons_bbtre_1000_2000 = [];  % HAB in [1000, 2000) m

profileCount = 0;

%% --- PROCESS BBTRE96 ---------------------------------------------------
for i = 1:75
    file96 = sprintf('%s%d.nc', base_path_96, i);
    [epsilons_bbtre_0_1000, epsilons_bbtre_1000_2000, used] = processFile_eps( ...
        file96, edges, epsilons_bbtre_0_1000, epsilons_bbtre_1000_2000, ...
        lon_min, lon_max, lat_min, lat_max, useRegionFilter);
    if used; profileCount = profileCount + 1; end 
end

%% --- PROCESS BBTRE97 ---------------------------------------------------
for i = 1:90
    file97 = sprintf('%s%d.nc', base_path_97, i);
    [epsilons_bbtre_0_1000, epsilons_bbtre_1000_2000, used] = processFile_eps( ...
        file97, edges, epsilons_bbtre_0_1000, epsilons_bbtre_1000_2000, ...
        lon_min, lon_max, lat_min, lat_max, useRegionFilter);
    if used; profileCount = profileCount + 1; end
end

%% --- SAVE --------------------------------------------------------------
save('eps_bbtre.mat', ...
     'epsilons_bbtre_0_1000', 'epsilons_bbtre_1000_2000');

fprintf('Done!   (%d profiles kept)\n', profileCount);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eps_0_1000, eps_1000_2000, usedProfile] = ...
    processFile_eps(filename, edges, ...
                    eps_0_1000, eps_1000_2000, ...
                    lon_min, lon_max, lat_min, lat_max, useRegionFilter)

usedProfile = false;

if ~isfile(filename)
    fprintf('File not found: %s\n', filename);
    return;
end

try
    % ---------- read data ----------------------------------------------
    lon      = ncread(filename, 'LONGITUDE');
    lat      = ncread(filename, 'LATITUDE');
    depth    = ncread(filename, 'DEPTH');      % [m]
    epsilon  = ncread(filename, 'EPSILON');    % [W/kg]
    botDepth = ncread(filename, 'BOT_DEPTH');  % [m]

    % ---------- if there is no record of bottom depth, skip-------------
    if isnan(botDepth) || botDepth <= 1000
        fprintf('Skipping %s: BOT_DEPTH = %.1f\n', filename, botDepth);
        return;
    end


    % ---------- if outside the region we care about, skip-------------
    if useRegionFilter
        if ~(lon >= lon_min && lon <= lon_max && ...
             lat >= lat_min && lat <= lat_max)
            % out of region
            return;
        end
    end

    % ---------- HAB & binning for epsilon ------------------------------
    HAB_eps = botDepth - depth;   % height above bottom for epsilon

    % bin 1: 0–1000 m HAB
    dmin1 = edges(1); dmax1 = edges(2);
    idx1  = HAB_eps >= dmin1 & HAB_eps < dmax1;
    if any(idx1)
        eps_0_1000 = [eps_0_1000; epsilon(idx1)];
    end

    % bin 2: 1000–2000 m HAB
    dmin2 = edges(2); dmax2 = edges(3);
    idx2  = HAB_eps >= dmin2 & HAB_eps < dmax2;
    if any(idx2)
        eps_1000_2000 = [eps_1000_2000; epsilon(idx2)];
    end

    usedProfile = true;

catch ME
    fprintf('Error reading %s → %s\n', filename, ME.message);
end
end
