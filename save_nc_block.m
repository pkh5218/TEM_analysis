function save_nc_block(filename, lat, lev, N_day, var_list)
    % SAVE_NC_BLOCK Save multiple 3D variables into NetCDF
    %   save_nc_block(filename, lat, lev, N_day, var_list)
    %   filename: output .nc filename
    %   lat, lev: latitude and pressure arrays
    %   N_day: number of time steps
    %   var_list: cell array, each row is {var_name, var_data}

    if isfile(filename)
        delete(filename);  % Avoid overwrite error
    end

    % Define dimensions
    nccreate(filename, 'lat', 'Dimensions', {'lat', length(lat)});
    nccreate(filename, 'lev', 'Dimensions', {'lev', length(lev)});
    nccreate(filename, 'day', 'Dimensions', {'day', N_day});
    ncwrite(filename, 'lat', lat);
    ncwrite(filename, 'lev', lev);
    ncwrite(filename, 'day', 1:N_day);
    ncwriteatt(filename, 'lat', 'units', 'degrees_north');
    ncwriteatt(filename, 'lev', 'units', 'hPa');
    ncwriteatt(filename, 'day', 'units', 'days since Jan 1');

    % Loop through all variables
    for i = 1:size(var_list,1)
        var_name = var_list{i,1};
        var_data = var_list{i,2};
        nccreate(filename, var_name, ...
            'Dimensions', {'lat', size(var_data,1), 'lev', size(var_data,2), 'day', N_day}, ...
            'DeflateLevel', 4);  % optional: compression
        ncwrite(filename, var_name, var_data);
    end
end

