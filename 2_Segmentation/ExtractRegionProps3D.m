%%
% Primordia Segmentation.
% Copyright (C) 2021 D. Eschweiler, J. Stegmaier
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the Liceense at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
% Please refer to the documentation for more information about the software
% as well as for installation instructions.
%
% If you use this application for your work, please cite the repository and one
% of the following publications:
%
% TBA
%
%%

addpath('../ThirdParty')
addpath('../ThirdParty/saveastiff_4.3');


%% specify the output directory
outputRoot = uigetdir('C:\Users\stegmaier\Downloads\GlobalOutputTest\', 'Please select the output root directory ...');
outputRoot = [outputRoot filesep];
if (~isfolder(outputRoot)); mkdir(outputRoot); end

%% specify the input directory
inputFolderTrackedNuclei = [outputRoot 'Nuclei_Tracked' filesep];
if (~isfolder(inputFolderTrackedNuclei))
    inputFolderTrackedNuclei = uigetdir('X:\Projects\KnautNYU_PrimordiumCellSegmentation\2021_10_08\20210506_Stiching_EGFP_Caax-H2a-mCherry\Nuclei_Tracked\', 'Please select the input directory (should contain tracked 3D tiffs for each frame)');
    inputFolderTrackedNuclei = [inputFolderTrackedNuclei filesep];
else
   disp(['Using tracked nuclei located in ' inputFolderTrackedNuclei]); 
end
inputFiles = dir([inputFolderTrackedNuclei '*.tif']);

%% specify the output directory
outputFolderRegionProps = [outputRoot 'Nuclei_Tracked_RegionProps' filesep];
if (~isfolder(outputFolderRegionProps)); mkdir(outputFolderRegionProps); end

zspacing = 2.4615; %% 0.325 x 0.325 x 0.8 µm 

%% setup the specifiers
specifiers = 'id;volume;xpos;ypos;zpos;equivDiameter;princialAxisLength1;principalAxisLength2;principalAxisLength3;orientation1;orientation2;orientation3;convexVolume;solidity;surfaceArea';

%% process all input files in parallel
parfor i=1:length(inputFiles)
    
    %% load the current image
    currentImage = loadtiff([inputFolderTrackedNuclei inputFiles(i).name]);
    
    %% resize the input image to get an isotropic one
    currentImage = imresize3(currentImage, [size(currentImage,1), size(currentImage,2), round(zspacing * size(currentImage,3))], 'nearest');
    
    %% extract the region props
    currentRegionProps = regionprops3(currentImage, 'Centroid', 'PrincipalAxisLength', 'Volume', 'SurfaceArea', 'ConvexVolume', 'EquivDiameter', 'Solidity', 'Solidity', 'Orientation');
    
    %% initialize the results table
    resultTable = table2array(currentRegionProps);
    
    %% specify feature indices
    volumeIndex = 1;
    centroidIndex = 2:4;
    equivDiameterIndex = 5;
    principalAxisLength = 6:8;
    orientationIndex = 9:11;
    convexVolumeIndex = 12;
    solidityIndex = 13;
    surfaceAreaIndex = 14;
    
    %% set the ids
    ids = 1:size(resultTable,1);
    resultTable = [ids', resultTable];
    resultTable(resultTable(:,2) == 0, :) = [];
    
    %% write the results to disk
    outputFileName = [outputFolderRegionProps strrep(inputFiles(i).name, '.tif', '_RegionProps.csv')];
    dlmwrite(outputFileName, resultTable, ';');
    prepend2file(specifiers, outputFileName, true);
end