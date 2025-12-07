clear; clc;


natre_path = '../../../../Data/NATRE/NATRE_data/Data_UCSD/NATRE/natre_';



dmin = 500;
dmax = 2000;



epsilons_natre_500_2000 = 0;

% --- PROCESS NATRE ---
for i = 1:155
    filename = sprintf('%s%d.nc', natre_path, i);
    if ~isfile(filename)
        fprintf('File not found: %s\n', filename);
        continue;
    end

    try
        lon      = ncread(filename, 'LONGITUDE');
        lat      = ncread(filename, 'LATITUDE');
        
        % Read the depth and epsilon data
        depth    = ncread(filename, 'DEPTH');    % 1D array of depths
        epsilon  = ncread(filename, 'EPSILON');  % 1D array of epsilon values
        

        % Bin epsilon using the original depth array
        inDepthRange = (depth >= dmin) & (depth < dmax);
        if any(inDepthRange)
            epsilons_natre_500_2000 = [epsilons_natre_500_2000; epsilon(inDepthRange)];
        end
       
    catch ME
        fprintf('Error reading file: %s\n  --> %s\n', filename, ME.message);
    end
end

% --- SAVE both the binned epsilon and N2 arrays ---
save('eps_natre.mat', 'epsilons_natre_500_2000');

