%% SG_to_netcdf
% Dump an SG_grid object to a NetCDF file
%
% Tristan Abbott // Massachusetts Institute of Technology // 08/21/2018
%
%%% Syntax
%   SG_to_netcdf(grid, fname)
%
%%% Description
% Dumps variables from an SG_grid object to a NetCDF file. 
% Note: currently implemented only for 3D scalar fields.
%
%%% Input Arguments
% *grid - model grid and fields:*
% grid must be created by SG_grid (or at least must have the format required for
% this function to run.
%
% *fileloc - path to the NetCDF file:*
%
function SG_to_netcdf(grid, fname)

    try
        ncid = netcdf.create(fname, 'NOCLOBBER');
    catch
        error('SG_to_netcdf:file_exists', ...
                sprintf('File %s already exists', fname));
    end

    % Create dimensions
    x_id = netcdf.defDim(ncid, 'x', grid.nx);
    y_id = netcdf.defDim(ncid, 'y', grid.ny);
    z_id = netcdf.defDim(ncid, 'z', grid.nzm);
    p_id = netcdf.defDim(ncid, 'p', grid.nzm);
    t_id = netcdf.defDim(ncid, 'time', 'NC_UNLIMITED');

    % Create variables
    fields = fieldnames(grid);
    var_id = [];
    var_name = {};
    for ii = 1:numel(fields)
        fld = grid.(fields{ii});
        if (size(fld,1) == grid.nzm & ...
           size(fld,2) == grid.ny & ...
           size(fld,3) == grid.nx)
            
            var_id = [var_id, ...
                netcdf.defVar(ncid, upper(fields{ii}), 'NC_DOUBLE', ...
                [x_id, y_id, z_id, t_id])];
            var_name{end+1} = fields{ii};
        end
    end

    netcdf.endDef(ncid);

    % Add data to dimensions
    netcdf.putVar(ncid, x_id, grid.x);
    netcdf.putVar(ncid, y_id, grid.y);
    netcdf.putVar(ncid, z_id, grid.z);
    netcdf.putVar(ncid, p_id, grid.p);

    % Add data to variables
    fields = fieldnames(grid);
    for ii = 1:numel(var_id)
        netcdf.putVar(ncid, var_id[ii], ...
            permute(grid.(var_name{ii}), [3 2 1]));
    end

    netcdf.close(ncid);

end
