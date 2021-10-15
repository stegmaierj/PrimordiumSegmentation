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

%% setup input paths
addpath('../ThirdParty')
addpath('../ThirdParty/saveastiff_4.3');

%% specify the output directory
outputFolder = uigetdir('C:\Users\stegmaier\Downloads\TGMM\', 'Please select a directory to store the results to ...');
outputFolder = [outputFolder filesep];
if (~isfolder(outputFolder)); mkdir(outputFolder); end

%% specify the input directory
inputFolder = uigetdir('X:\Projects\KnautNYU_PrimordiumCellSegmentation\2021_10_08\20210506_Stiching_EGFP_Caax-H2a-mCherry\Nuclei_Tracked\', 'Please select the input directory (should contain tracked 3D tiffs for each frame)');
inputFolder = [inputFolder filesep];
inputFiles = dir([inputFolder '*.tif']);

%% define z-spacing (only needed if isotropic images should be processed).
zspacing = 1; %% 0.325 x 0.325 x 0.8 µm 

%% specify indices of the features in the region props table
volumeIndex = 1;
centroidIndex = 2:4;
equivDiameterIndex = 5;
principalAxisLength = 6:8;
orientationIndex = 9:11;
convexVolumeIndex = 12;
solidityIndex = 13;
surfaceAreaIndex = 14;

%% initialize the previous labels array
previousLabels = [];

%% process all image files sqeuentially
for i=1:length(inputFiles)
    
    %% read the current input image
    currentImage = loadtiff([inputFolder inputFiles(i).name]);
    
    %% potentially resize the input image in z if needed
    if (zspacing ~= 1)
        currentImage = imresize3(currentImage, [size(currentImage,1), size(currentImage,2), round(zspacing * size(currentImage,3))], 'nearest');
    end
    
    %% extract the region properties
    currentRegionProps = regionprops3(currentImage, 'Centroid', 'PrincipalAxisLength', 'Volume', 'SurfaceArea', 'ConvexVolume', 'EquivDiameter', 'Solidity', 'Solidity', 'Orientation');
    
    %% convert region props table to a regular array
    resultTable = table2array(currentRegionProps);
    
    %% create ids for the current list of objects
    ids = 1:size(resultTable,1);
    
    %% assemble the results table
    resultTable = [ids', resultTable(:,volumeIndex), round(resultTable(:,centroidIndex)), (i-1)*ones(size(ids')), 100*ones(size(ids')), ids', ids'];
    resultTable(resultTable(:,2) == 0, :) = [];
        
    %% create TGMM style output xml files for each of the frames
    docNode = com.mathworks.xml.XMLUtils.createDocument('document'); %#ok<JAPIMATHWORKS>
    toc = docNode.getDocumentElement;
    
    %% add each detection as a fake Gaussian mixture model
    %% only the centroid and IDs are reasonable values, 
    %% the rest are just dummy values from one of the TGMM examples.
    for j=1:size(resultTable,1)
        currentGM = docNode.createElement('GaussianMixtureModel');
        currentGM.setAttribute('id', num2str(resultTable(j,1)));
        currentGM.setAttribute('lineage', num2str(resultTable(j,1)));
        
        if (ismember(resultTable(j,1), previousLabels))
            currentGM.setAttribute('parent', num2str(resultTable(j,1)));
        else
            currentGM.setAttribute('parent', '-1');
        end
        currentGM.setAttribute('dims', '3');
        currentGM.setAttribute('splitScore', '3');
        currentGM.setAttribute('scale', '1 1 1');
        currentGM.setAttribute('nu', '126.937');
        currentGM.setAttribute('beta', '126.937');
        currentGM.setAttribute('m', num2str(resultTable(j,3:5)));
        currentGM.setAttribute('W', '1 0 0 0 1 0 0 0 1');
        currentGM.setAttribute('nuPrior', '4');
        currentGM.setAttribute('betaPrior', '0.075017');
        currentGM.setAttribute('alphaPrior', '0');
        currentGM.setAttribute('distMRFPrior', '0');
        currentGM.setAttribute('mPrior', num2str(resultTable(j,3:5)));
        currentGM.setAttribute('WPrior', '1 0 0 0 1 0 0 0 1');
        currentGM.setAttribute('svIdx', '0');
        toc.appendChild(currentGM);
    end

    %% write the current result file
    outputFileName = [outputFolder strrep(inputFiles(i).name, '.tif', '_RegionProps.xml')];
    xmlwrite(outputFileName,docNode);
    
    %% store the previous labels to see if there's a predecessor available    
    previousLabels = resultTable(:,1);
end