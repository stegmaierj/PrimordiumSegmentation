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

% %% setup input paths
% addpath('../ThirdParty')
% addpath('../ThirdParty/saveastiff_4.3');
% addpath('../ThirdParty/matlab-mastodon-importer/src');
% 
% %% specify the output directory
% outputRoot = uigetdir('C:\Users\stegmaier\Downloads\GlobalOutputTest\', 'Please select the output root directory ...');
% outputRoot = [outputRoot filesep];
% if (~isfolder(outputRoot)); mkdir(outputRoot); end
% 
% %% specify the output directory
% outputFolderCorrected = [outputRoot 'Nuclei_Tracked_Corrected' filesep];
% if (~isfolder(outputFolderCorrected)); mkdir(outputFolderCorrected); end
% 
% %% create subfolders for csv and image output
% outputFolderImages = [outputFolderCorrected 'TrackedImages' filesep];
% outputFolderCSV = [outputFolderCorrected 'RegionProps' filesep];
% if (~isfolder(outputFolderImages)); mkdir(outputFolderImages); end
% if (~isfolder(outputFolderCSV)); mkdir(outputFolderCSV); end
% 
% %% 
% inputFolderNuclei = [outputRoot 'Nuclei' filesep];
% if (~isfolder(inputFolderNuclei))
%     inputFolderNuclei = uigetdir('X:\Projects\KnautNYU_PrimordiumCellSegmentation\2021_10_08\20210506_Stiching_EGFP_Caax-H2a-mCherry\Nuclei\', 'Please select the input directory (should contain tracked 3D tiffs for each frame)');
%     inputFolderNuclei = [inputFolderNuclei filesep];
% else
%     disp(['Using raw nuclei image files contained in folder ' inputFolderNuclei]);
% end
% inputFilesRaw = dir([inputFolderNuclei '*.tif']);
% 
% %% specify the input directory
% inputFolderTrackedNuclei = [outputRoot 'Nuclei_Tracked' filesep];
% if (~isfolder(inputFolderTrackedNuclei))
%     inputFolderTrackedNuclei = uigetdir('X:\Projects\KnautNYU_PrimordiumCellSegmentation\2021_10_08\20210506_Stiching_EGFP_Caax-H2a-mCherry\Nuclei_Tracked\', 'Please select the input directory (should contain tracked 3D tiffs for each frame)');
%     inputFolderTrackedNuclei = [inputFolderTrackedNuclei filesep];
% else
%     disp(['Using tracked image files contained in folder ' inputFolderTrackedNuclei]);
% end
% inputFilesSegmentation = dir([inputFolderTrackedNuclei '*.tif']);
% 
% %% specify the input mastodon file
% [inputFileMastodon, inputPathMastodon]  = uigetfile('*.mastodon', 'Select the mastodon file you want to apply to the selected images');
% 
% %% bring mastodon file into a d_orgs like format for easier application to the image data
% [G, metadata, tss] = import_mastodon( [inputPathMastodon inputFileMastodon] );
% 
% nodes = table2array(G.Nodes(:,1:5));
% nodes(:,1) = nodes(:,1)+1;
% nodes(:,5) = nodes(:,5)+1;
% 
% edges = table2array(G.Edges);
% 
% %% sort the edges based on their start time point
% %% otherwise it can happen that objects are added later and start a new track
% edges(:,end+1) = 0;
% for i=1:size(edges,1)
%     
%     currentEdge = edges(i,:);
%     
%     object1 = nodes(currentEdge(1), :);
%    
%     timePoint1 = object1(5);
%     edges(i,end) = timePoint1;
% end
% 
% [timePoints, edgeOrder] = sort(edges(:,end), 'ascend');
% edges = edges(edgeOrder, :);
% 
% visitedNodes = zeros(size(nodes,1),1);
% visitedEdges = zeros(size(edges,1),1);
% 
% d_orgs = zeros(size(nodes,1), 6);
% for i=1:size(nodes,1)
%     currentTimePoint = nodes(i,5);
%     d_orgs(i, currentTimePoint, 1) = nodes(i,1);
%     d_orgs(i, currentTimePoint, 2) = 1;
%     d_orgs(i, currentTimePoint, 3:5) = round(nodes(i,2:4));
% end
% 
% deletionIndices = zeros(size(d_orgs,1),1);
% for i=1:size(edges,1)
%     currentEdge = edges(i,:);
%     
%     object1 = nodes(currentEdge(1), :);
%     object2 = nodes(currentEdge(2), :);
%     
%     timePoint1 = object1(5);
%     timePoint2 = object2(5);
%     
%     object1Index = find(d_orgs(:,timePoint1,1) == currentEdge(1));
%     object2Index = find(d_orgs(:,timePoint2,1) == currentEdge(2));
%         
%     d_orgs(object1Index, timePoint2, :) = d_orgs(object2Index, timePoint2, :);
%     d_orgs(object2Index, timePoint2, :) = 0;
%     
%     deletionIndices(object2Index) = 1;
% end
% 
% %% remove empty rows
% d_orgs(deletionIndices > 0, :, :) = [];
% 
% %% renumber indices that belong to a single track
% for i=1:size(d_orgs,1)
%     validIndices = squeeze(d_orgs(i,:,1)) > 0;
%     d_orgs(i,validIndices,1) = i;
% end
% 
% %% define z-spacing (only needed if isotropic images should be processed).
% zspacing = 1; %% 0.325 x 0.325 x 0.8 µm
% 
% %% specify indices of the features in the region props table
% volumeIndex = 1;
% centroidIndex = 2:4;
% equivDiameterIndex = 5;
% principalAxisLength = 6:8;
% orientationIndex = 9:11;
% convexVolumeIndex = 12;
% solidityIndex = 13;
% surfaceAreaIndex = 14;
% 
% %% initialize the previous labels array
% previousLabels = [];
% 
% %% save options
% clear options;
% options.compress = 'lzw';
% options.overwrite = true;

%% process all image files sqeuentially
for i=338%1:length(inputFilesSegmentation)
    
    %% read the current input image
    currentSegmentationImage = loadtiff([inputFolderTrackedNuclei inputFilesSegmentation(i).name]);
    
    %% read the current raw image
    currentRawImage = loadtiff([inputFolderNuclei inputFilesRaw(i).name]);
    
    %% potentially resize the input image in z if needed
    if (zspacing ~= 1)
        currentSegmentationImage = imresize3(currentSegmentationImage, [size(currentSegmentationImage,1), size(currentSegmentationImage,2), round(zspacing * size(currentSegmentationImage,3))], 'nearest');
    end
    
    %% create a seed image
    currentSeedImage = zeros(size(currentSegmentationImage));
    validIndices = find(d_orgs(:,i,1) > 0);
    for j=validIndices'
        currentCentroid = round(squeeze(d_orgs(j,i,3:5)));
        currentSeedImage(currentCentroid(2), currentCentroid(1), currentCentroid(3)) = j;
    end
    
    %% extract the region properties
    currentRegionProps = regionprops3(currentSegmentationImage, 'Centroid', 'PrincipalAxisLength', 'Volume', 'SurfaceArea', 'ConvexVolume', 'EquivDiameter', 'Solidity', 'Solidity', 'Orientation', 'VoxelIdxList', 'BoundingBox');
    
    visitedIndices = zeros(size(currentRegionProps, 1), 1);
    
    resultImage = zeros(size(currentSegmentationImage));
    for j=validIndices'
        
        currentCentroid = round(squeeze(d_orgs(j,i,3:5)));
        oldLabel = currentSegmentationImage(currentCentroid(2), currentCentroid(1), currentCentroid(3));
                
        if (oldLabel == 0)
            resultImage(currentCentroid(2), currentCentroid(1), currentCentroid(3)) = j;
            continue;
        else
            visitedIndices(oldLabel) = 1;
        end
        
        %% get the labels of the current region
        currentLabels = unique(currentSeedImage(cell2mat(currentRegionProps(oldLabel, {'VoxelIdxList'}).VoxelIdxList)));
        numLabels = length(currentLabels);
        
        %% handle cases where a former segment contains more than one label
        if (numLabels > 2)
            
            borderPadding = 2;
            currAABB = round(currentRegionProps(oldLabel, {'BoundingBox'}).BoundingBox);
            rangeY = max(1, currAABB(1)-borderPadding):min(size(currentSegmentationImage,2), (currAABB(1)+currAABB(4)+borderPadding));
            rangeX = max(1, currAABB(2)-borderPadding):min(size(currentSegmentationImage,1), (currAABB(2)+currAABB(5)+borderPadding));
            rangeZ = max(1, currAABB(3)-borderPadding):min(size(currentSegmentationImage,3), (currAABB(3)+currAABB(6)+borderPadding));
            
            %% perform the split
            currentMask = double(currentSegmentationImage(rangeX, rangeY, rangeZ) == oldLabel);
            currentBoundary = imdilate(currentMask, strel('sphere', 1)) - currentMask;
            currentRegion = imgaussfilt3(double(currentRawImage(rangeX, rangeY, rangeZ)), 1) .* currentMask;
            currentRegion = max((max(currentRegion(:)) - currentRegion)  .* currentMask, max(currentRegion(:)) * currentBoundary);
            
            labelImage = imdilate(currentSeedImage(rangeX, rangeY, rangeZ) .* currentMask, strel('sphere', 1));
            currentRegion = uint16(imimposemin(currentRegion, labelImage > 0));
            
            currenRegionWs = watershed(currentRegion);
            splitImage = double(currenRegionWs > 0) .* double(currentSegmentationImage(rangeX, rangeY, rangeZ) == oldLabel);
            figure(2); imagesc(splitImage(:,:,round(size(splitImage,3)/2)));
            
            %% replace the original segment with the split one
            %% TODO
            oldRegionContent = currentSegmentationImage(rangeX, rangeY, rangeZ);
            oldRegionContent(oldRegionContent == oldLabel) = 0;
                        
            %% set tracking label for the split segment
            localRegionProps = regionprops(splitImage > 0, 'PixelIdxList');
            for k=1:length(localRegionProps)
                currentLabel = max(labelImage(localRegionProps(k).PixelIdxList));
                oldRegionContent(localRegionProps(k).PixelIdxList) = currentLabel;
            end
            
            %% as a fallback option, also add the labels (in case the watershed failed and extracted segments don't match seeds
            for k=currentLabels'
               if (k==0)
                   continue;
               else
                  oldRegionContent(labelImage == k) = k;
               end
            end
            
            resultImage(rangeX, rangeY, rangeZ) = oldRegionContent;
            
        else
            %% everything alright, just copy the segment with the new id
            resultImage(cell2mat(currentRegionProps(oldLabel, {'VoxelIdxList'}).VoxelIdxList)) = j;
        end
    end
    
    testStats = visitedIndices == (table2array(currentRegionProps(:,1)) > 0);
    
    errors = length(testStats) - sum(testStats)
    
    %% save the result image
    saveastiff(uint16(resultImage), [outputFolderImages inputFilesSegmentation(i).name], options);
end
