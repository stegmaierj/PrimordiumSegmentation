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

%% obtain input and output paths
flowSegDir = uigetdir(pwd, 'Select result folder containing the flow-based segmentation!');
wsSegDir = uigetdir(pwd, 'Select result folder containing the watershed-based segmentation!');
flowSegFiles = dir([flowSegDir filesep '*.tif']);
wsSegFiles = dir([wsSegDir filesep '*.tif']);

%% specify the output directory
outputRoot = uigetdir('C:\Users\stegmaier\Downloads\20210715_Stiching_EGFP_Caax-H2a-mCherry_crop\', 'Please select the output root directory ...');
outputFolder = [outputRoot filesep 'Nuclei_Segmented' filesep];
if (~isfolder(outputFolder)); mkdir(outputFolder); end

%% specify parameters
hMaxValue = 1;
minArea = 100;

%% get the number of images
numImages = length(flowSegFiles);

%% specify tiff writing options
clear options;
options.overwrite = true;
options.compress = 'lzw';

%% process all images in parallel
parfor f=1:numImages
    
    %% load the current input images    
    segImageFlow = loadtiff([flowSegDir filesep flowSegFiles(f).name]);
    segImageWs = loadtiff([wsSegDir filesep wsSegFiles(f).name]);
    
    %% obtain the image size and initialize result image
    imageSize = size(segImageFlow);
    resImage = zeros(imageSize);
    
    %% obtain region props
    regionProps = regionprops(segImageFlow, 'Area', 'BoundingBox', 'PixelIdxList');
    
    %% process all segments iteratively and optionally split them
    currentLabel = 1;
    for i=1:length(regionProps)
        
        %% skip processing small elements
        if (regionProps(i).Area <= minArea)
            continue;
        end
        
        %% convert aabb coordinates to index ranges        
        aabb = round(regionProps(i).BoundingBox);
        rangeX = max(1, aabb(1)):min((aabb(1)+aabb(4)), imageSize(2));
        rangeY = max(1, aabb(2)):min((aabb(2)+aabb(5)), imageSize(1));
        rangeZ = max(1, aabb(3)):min((aabb(3)+aabb(6)), imageSize(3));
        
        %% get image snippets of the current region
        flowSeg = segImageFlow(rangeY, rangeX, rangeZ) == i;
        watershedSeg = bwlabeln(segImageWs(rangeY, rangeX, rangeZ));

        %% identify existing labels
        validIndices = find(flowSeg > 0);
        labels = unique(watershedSeg(validIndices));

        %% check if splitting is required
        if (sum(labels > 0) > 1)
            
            %% create distance map
            currentDistMap = bwdist(~flowSeg);
            currentDistMap = max(currentDistMap(:)) - currentDistMap;
            
            %% create seeded segmentation input image and apply watershed
            hminImage = imimposemin(currentDistMap, (watershedSeg > 0) .* flowSeg);
            watershedImage = watershed(hminImage);
            
            %% apply watershed lines to the previous segmentation
            resultSnippet = flowSeg .* (watershedImage > 0);
            
            %% fill the current snippet with the split segment including new labels
            resultRegionProps = regionprops(resultSnippet > 0, 'PixelIdxList');
            currentImageRes = zeros(size(resultSnippet));
            for j=1:length(resultRegionProps)
                currentImageRes(resultRegionProps(j).PixelIdxList) = currentLabel;
                currentLabel = currentLabel + 1;
            end
            
            %% add segments to the result image
            resImage(rangeY, rangeX, rangeZ) = max(resImage(rangeY, rangeX, rangeZ), currentImageRes);
            
        %% if no split is required, just copy the previous segment
        else
            resImage(rangeY, rangeX, rangeZ) = max(resImage(rangeY, rangeX, rangeZ), flowSeg * currentLabel);
            currentLabel = currentLabel + 1;
        end
        
    end
        
    %% write the cleaned result image
    saveastiff(uint16(resImage), [outputFolder filesep flowSegFiles(f).name], options);
end