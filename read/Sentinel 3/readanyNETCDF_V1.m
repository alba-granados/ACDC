%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%    isardSAT S.L.   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% readanyNETCDF reads any file in NETCDF format
% INPUT : filename
% OUTPUT: a Matlab structure with all the data and attributes from the NETCDF file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alldata]=readanyNETCDF(file)

ncid = netcdf.open(file,'NC_NOWRITE'); %open file
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid); %get global attributes

for i=0:(nvars-1)
    [data(i+1).varname,data(i+1).xtype,data(i+1).dimids,data(i+1).natts] = netcdf.inqVar(ncid,i); %get variables name
    alldata.data.(data(i+1).varname) = netcdf.getVar(ncid,netcdf.inqVarID(ncid,data(i+1).varname)); %get variable data
    for j=0:(data(i+1).natts - 1)
        att.name = netcdf.inqAttName(ncid,netcdf.inqVarID(ncid,data(i+1).varname),j); %get attributes name
        att.name2put = att.name;
        if att.name(1) == '_'
           att.name2put(1) = []; 
        end
        alldata.attributes.(data(i+1).varname).(att.name2put) = netcdf.getAtt(ncid,netcdf.inqVarID(ncid,data(i+1).varname),(att.name)); %get attributes info
    end
end;

netcdf.close(ncid);

end

% [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid)            Return information about netCDF file

% [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,varid)      Return information about variable
% data = netcdf.getVar(ncid,varid)                              Return data from netCDF variable
% varid = netcdf.inqVarID(ncid,varname)                         Return ID associated with variable name

% attname = netcdf.inqAttName(ncid,varid,attnum)                Return name of netCDF attribute
% [xtype,attlen] = netcdf.inqAtt(ncid,varid,attname)            Return information about netCDF attribute
% attrvalue = netcdf.getAtt(ncid,varid,attname)                 Return netCDF attribute

