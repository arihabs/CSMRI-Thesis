function [imaVOL,scaninfo] = myloadminc(filename)
%function [imaVOL,scaninfo] = myloadminc(filename)
% 
% Modified version of loadminc.m
%==========================================================================
% Original Description of loadminc.m
%
% Function to load minc format input file. 
% This function use the netcdf MATLAB utility
%
% Matlab library function for MIA_gui utility. 
% University of Debrecen, PET Center/LB 2010
%==========================================================================
% Modified and renamed by Ariel Habshush.
% 
% Description:
% This function makes changes to the loadminc.m function. The original
% function scaled/normalized the image/volume data so that pixel
% values no longer had unsigned integer bit values. This scaling has
% been removed.
%--------------------------------------------------------------------------
if nargin == 0
     [FileName, FilePath] = uigetfile('*.mnc','Select minc file');
     filename = [FilePath,FileName];
     if FileName == 0;
          imaVOL = [];scaninfo = [];
          return;
     end
end

ncid=netcdf.open(filename,'NC_NOWRITE');
scaninfo.filename = filename;

% Get voxel size.
pixsizex = netcdf.getAtt(ncid,netcdf.inqVarID(ncid,'xspace'),'step');
pixsizey = netcdf.getAtt(ncid,netcdf.inqVarID(ncid,'yspace'),'step');
pixsizez = netcdf.getAtt(ncid,netcdf.inqVarID(ncid,'zspace'),'step');
scaninfo.pixsize = abs([pixsizex pixsizey pixsizez]);

% Get units of voxel size.
unitX = netcdf.getAtt(ncid,netcdf.inqVarID(ncid,'xspace'),'units');
unitY = netcdf.getAtt(ncid,netcdf.inqVarID(ncid,'yspace'),'units');
unitZ = netcdf.getAtt(ncid,netcdf.inqVarID(ncid,'zspace'),'units');
scaninfo.pixunits = [unitX, unitY, unitZ];

% Get starting location of scan.
x_start = netcdf.getAtt(ncid,netcdf.inqVarID(ncid,'xspace'),'start');
if isempty(x_start)
    x_start = 0;
end
y_start = netcdf.getAtt(ncid,netcdf.inqVarID(ncid,'yspace'),'start');
if isempty(y_start)
    y_start = 0;
end
z_start = netcdf.getAtt(ncid,netcdf.inqVarID(ncid,'zspace'),'start');
if isempty(z_start)
    z_start = 0;
end

scaninfo.space_start = ([x_start y_start z_start]);

% Get actual volume data.
varid = netcdf.inqVarID(ncid,'image');
volume = netcdf.getVar(ncid,varid);

% When reading in image data, the y-dimension is the first dimension,
% and the x-dimension is the 2nd dimension of the volume array. These
% dimesions are now changed so that the x-dimension is the 1st
% dimension and the y-dimension is the 2nd dimension.

imaVOL = permute(volume,[2,1,3]);

% Get number of slices.
scaninfo.num_of_slice  = size(imaVOL,3);

% Get Image Dimension.
scaninfo.ImageDim = [size(imaVOL,1) size(imaVOL,2)];
scaninfo.Frames = 1;
scaninfo.start_times = [];
scaninfo.tissue_ts = [];
scaninfo.frame_lengths = [];
scaninfo.FileType    = 'mnc';

% Close mnc file.
netcdf.close(ncid);

end %function



