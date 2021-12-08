% mnc2mat.m
% 
% This script loads all .mnc files (containing 3D volume data) in a
% specified location into 3D arrays and stores the array as a .mat
% file.
%--------------------------------------------------------------------------
addpath('../MATLABpackages/');
dataLocation = 'F:\Thesis\Thesis_Datasets\mcguill\Brainweb\TwentyNormalBrains';
subfolders = dir(dataLocation);
% Filter out all the items in the main folder that are not
% directories.
subfolders(~[subfolders.isdir]) = [];
% Filters out the parent and current directory '.' and '..'
tf = ismember( {subfolders.name}, {'.', '..'});
subfolders(tf) = [];
numberOfFolders = length(subfolders);

% Loop through Subfolders
for iii = 1:numberOfFolders
    subfolderPath = fullfile(dataLocation, subfolders(iii).name);
    mncfile = dir(fullfile(subfolderPath, '*.mnc'));
    fullFileName = fullfile( subfolderPath, mncfile.name );
    [imageVolume, scaninfo] = myloadminc(fullFileName);
    
    % Save the data with the same name as the mnc file without the
    % file extension.
    [pathstr, name, ext] = fileparts(fullFileName);
    saveFileName = fullfile(pathstr, name);
    
    save(saveFileName, 'imageVolume', 'scaninfo');
    
    
end
