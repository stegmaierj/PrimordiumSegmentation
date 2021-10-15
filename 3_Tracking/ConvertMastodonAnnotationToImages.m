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
addpath('../ThirdParty/matlab-mastodon-importer/src');

%% specify the output directory
outputFolder = uigetdir('C:\Users\stegmaier\Downloads\MastodonImages\', 'Please select a directory to store the results to ...');
outputFolder = [outputFolder filesep];
if (~isfolder(outputFolder)); mkdir(outputFolder); end

%% create subfolders for csv and image output
outputFolderImages = [outputFolder 'TrackedImages' filesep];
outputFolderCSV = [outputFolder 'RegionProps' filesep];

%% specify the input directory
inputFolder = uigetdir('X:\Projects\KnautNYU_PrimordiumCellSegmentation\2021_10_08\20210506_Stiching_EGFP_Caax-H2a-mCherry\Nuclei_Tracked\', 'Please select the input directory (should contain tracked 3D tiffs for each frame)');
inputFolder = [inputFolder filesep];
inputFiles = dir([inputFolder '*.tif']);

%% specify the input mastodon file
[inputFileMastodon, inputPathMastodon]  = uigetfile('*.mastodon', 'Select the mastodon file you want to apply to the selected images');

%% bring mastodon file into a d_orgs like format for easier application to the image data
[ G, metadata, tss ] = import_mastodon( [inputPathMastodon inputFileMastodon] );

nodes = table2array(G.Nodes(:,1:5));
nodes(:,1) = nodes(:,1)+1;
nodes(:,5) = nodes(:,5)+1;

edges = table2array(G.Edges);

visitedNodes = zeros(size(nodes,1),1);
visitedEdges = zeros(size(edges,1),1);

d_orgs = zeros(size(nodes,1), 6);
for i=1:size(nodes,1)
   currentTimePoint = nodes(i,5);   
   d_orgs(i, currentTimePoint, 1) = nodes(i,1);
   d_orgs(i, currentTimePoint, 2) = 1;
   d_orgs(i, currentTimePoint, 3:5) = nodes(i,2:4);     
end

deletionIndices = zeros(size(d_orgs,1),1);
for i=1:size(edges,1)
   currentEdge = edges(i,:);
   
   object1 = nodes(currentEdge(1), :);
   object2 = nodes(currentEdge(2), :);
   
   timePoint1 = object1(5);
   timePoint2 = object2(5);
   
   object1Index = find(d_orgs(:,timePoint1,1) == currentEdge(1));
   object2Index = find(d_orgs(:,timePoint2,1) == currentEdge(2));
   
   d_orgs(object1Index, timePoint2, :) = d_orgs(object2Index, timePoint2, :);
   d_orgs(object2Index, timePoint2, :) = 0;
   
   deletionIndices(object2Index) = 1;
end

%% remove empty rows
d_orgs(deletionIndices > 0, :, :) = [];

%% renumber indices that belong to a single track
for i=1:size(d_orgs,1)
   validIndices = squeeze(d_orgs(i,:,1)) > 0;
   d_orgs(i,validIndices,1) = i;
end

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

%% save options
clear options;
options.compress = 'lzw';
options.overwrite = true;

%% process all image files sqeuentially
for i=1:length(inputFiles)
    
    %% read the current input image
    currentImage = loadtiff([inputFolder inputFiles(i).name]);
    
    %% potentially resize the input image in z if needed
    if (zspacing ~= 1)
        currentImage = imresize3(currentImage, [size(currentImage,1), size(currentImage,2), round(zspacing * size(currentImage,3))], 'nearest');
    end
    
    %% extract the region properties
    currentRegionProps = regionprops3(currentImage, 'Centroid', 'PrincipalAxisLength', 'Volume', 'SurfaceArea', 'ConvexVolume', 'EquivDiameter', 'Solidity', 'Solidity', 'Orientation', 'VoxelIdxList');
    
    validIndices = find(d_orgs(:,i,1) > 0);
    resultImage = zeros(size(currentImage));
    for j=validIndices'
        
        currentCentroid = squeeze(d_orgs(j,i,3:5));
        oldLabel = currentImage(currentCentroid(2), currentCentroid(1), currentCentroid(3));
        
        if (oldLabel == 0)
            disp('Potential error detected!!');
            
            %% TODO: how to handle cells where the centroid is outside of the cell?? Obviously an undersegmentation error.
%                         
%             figure(1); clf;
%             imagesc(currentImage(:,:,currentCentroid(3))); hold on;
%             plot(currentCentroid(1), currentCentroid(2), '*r');

            continue;
        end
        
        
        resultImage(cell2mat(currentRegionProps(oldLabel, {'VoxelIdxList'}).VoxelIdxList)) = j;
        
        %% TODO: handle splits!
    end
    
    %% save the result image
    saveastiff(uint16(resultImage), [outputFolderImages inputFiles(i).name], options);
end
