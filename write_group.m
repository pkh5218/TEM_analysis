function write_group(outfile, lat, lev_hPa, NDays, name_data_pairs)
    if exist(outfile,'file'), delete(outfile); end

    nccreate(outfile,'lat', 'Dimensions',{'lat',numel(lat)}, 'Datatype','double');
    ncwrite (outfile,'lat', lat);  ncwriteatt(outfile,'lat','units','degrees_north');

    nccreate(outfile,'lev', 'Dimensions',{'lev',numel(lev_hPa)}, 'Datatype','double');
    ncwrite (outfile,'lev', lev_hPa); ncwriteatt(outfile,'lev','units','hPa');

    nccreate(outfile,'day', 'Dimensions',{'day',NDays}, 'Datatype','double');
    ncwrite (outfile,'day', 1:NDays);   ncwriteatt(outfile,'day','units','days since Jan-1');

    for k = 1:size(name_data_pairs,1)
        vname = name_data_pairs{k,1};
        vdata = name_data_pairs{k,2};
        nccreate(outfile, vname, 'Dimensions', {'lat', numel(lat), 'lev', numel(lev_hPa), 'day', NDays}, ...
                 'DeflateLevel', 3, 'Shuffle', true, 'Format','netcdf4');
        ncwrite(outfile, vname, vdata(:,:,1:NDays));
    end
end