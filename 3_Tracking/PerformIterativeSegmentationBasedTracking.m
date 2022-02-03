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

%% Performs segmentation based on registration and clustering.

%% load additional dependencies
addpath('../ThirdParty/saveastiff_4.3/');

%%%%%%%%%% PARAMETERS %%%%%%%%%%%
useRegistrationBasedCorrection = false;
preloadData = true; %% only set to false if debugging to avoid reloading the data at each call
minVolume = 100; %% excludes objects smaller than this volume
%%%%%%%%%% PARAMETERS %%%%%%%%%%%

%% specify elastix path
elastixRootDir = [pwd '/../ThirdParty/Elastix/'];
if (ismac)
    elastixPath = [elastixRootDir 'macOS/bin/'];
elseif (ispc)
    elastixPath = [elastixRootDir 'Windows/'];
    elastixPath = strrep(elastixPath, '/', filesep);
else
    elastixPath = [elastixRootDir 'Linux/bin/'];
end

%% specify the output directory
outputRoot = uigetdir('C:\Users\stegmaier\Downloads\20211129_Stiching_EGFP_Caax-H2a-mCherry_CA-Mypt1\', 'Please select the output root directory ...');
outputRoot = [outputRoot filesep];
if (~isfolder(outputRoot)); mkdir(outputRoot); end

%% get the input directories
inputDir = [outputRoot 'Nuclei_Segmented' filesep];
if (~isfolder(inputDir))
    inputDir = uigetdir('I:\Projects\2021\KnautNYU_PrimordiumCellSegmentation\Processing\item_0012_GradientVectorFlowTrackingImageFilter\', 'Specify the input folder containing segmentation tiffs for each time point.');
    inputDir = [inputDir filesep];
else
   disp(['Using segmented images contained in the folder ' inputDir]);    
end

%% get the transformations directory
if (useRegistrationBasedCorrection == true)
    transformationDir = [outputRoot 'Transformations' filesep];
    if (~isfolder(transformationDir))
        transformationDir = uigetdir('D:\ScieboDrive\Projects\2021\KnautNYU_PrimordiumCellSegmentation\Processing\Transformations\', 'Specify the transformation directory');
        transformationDir = [transformationDir filesep];
        transformationFiles = dir([transformationDir '*.txt']);
    else
        disp(['Using transformations contained in the folder ' transformationDir]);    
    end
end

%% create the result folder for the tracking data
outputFolderTracked = [outputRoot 'Nuclei_Tracked' filesep];
if (~isfolder(outputFolderTracked)); mkdir(outputFolderTracked); end

%% get information about the input/output files
inputFiles = dir([inputDir '*.tif']);
numFrames = length(inputFiles);

%% specify save options
clear options;
options.overwrite = true;
options.compress = 'lzw';

%% preload data to have all frames in memory for faster processing
if (preloadData == true)
    
    %% clear previous data
    clear resultImages;
    clear rawImages;
    clear segmentationImages;
    clear regionProps;
    
    %% load and transform all input images
    parfor i=1:numFrames
        
        %% load the current image and obtain the region props
        segmentationImages{i} = loadtiff([inputDir inputFiles(i).name]);
        regionProps{i} = regionprops(segmentationImages{i}, 'Area', 'Centroid', 'PixelIdxList', 'BoundingBox');
                
        %% remove spurious objects that have their centroid outside of the valid area
        for j=1:size(regionProps{i})
            if (regionProps{i}(j).Area <= 0)
                continue;
            end
            
            currentCentroid = round(regionProps{i}(j).Centroid);
            currentLabel = segmentationImages{i}(currentCentroid(2), currentCentroid(1), currentCentroid(3));
            
            if (currentLabel ~= j)
                segmentationImages{i}(regionProps{i}(j).PixelIdxList) = 0;
            end
        end
        
        %% recompute the region props after removing spurious objects
        regionProps{i} = regionprops(segmentationImages{i}, 'Area', 'Centroid', 'PixelIdxList', 'BoundingBox');
        
        %% transform the image at time point t to time point t-1 for better overlaps
        if (useRegistrationBasedCorrection == true)
            
            %% save temporary image
            outputDirTemp = [tempdir num2str(i) filesep];
            inputFile = [outputDirTemp 'currImage.tif'];
            saveastiff(uint16(segmentationImages{i}), inputFile, options);
            
            %% apply transformation to the segmentation image
            transformixCommand = [elastixPath 'transformix.sh ' ...
                '-in ' inputFile ' ' ...
                '-out ' outputDirTemp ' ' ...
                '-tp ' transformationDir transformationFiles(i).name];
            
            if (ispc)
                transformixCommand = strrep(transformixCommand, 'transformix.sh', 'transformix.exe');
            end
            system(transformixCommand);
            
            %% load the transformed image at time point t and the previous image at t-1
            transformedCurrentImage = loadtiff([outputDirTemp 'result.tif']);
            movefile([outputDirTemp 'result.tif'], [resultDirTrans strrep(inputFiles(i).name, '.tif', '_Trans.tif')]); rmdir(outputDirTemp, 's');
            
            %% save the regionprops for later computations
            regionPropsTransformed{i} = regionprops(transformedCurrentImage, 'Area', 'Centroid', 'PixelIdxList', 'BoundingBox');
        else
            regionPropsTransformed{i} = regionProps{i};
        end
        
        %% initialize the result images
        resultImages{i} = uint16(zeros(size(segmentationImages{i})));
    end
end

%% clear previous tracking results
clear visitedIndices;
for i=1:numFrames
    visitedIndices{i} = zeros(size(regionProps{i},1)); %#ok<SAGROW>
    resultImages{i}(:) = 0;
end

%% initialize current tracking id
currentTrackingId = 1;

%% iterate over all frames in a backward fashion
for i=numFrames:-1:1
    
    %% process all cells contained in the current frame
    for j=1:size(regionProps{i},1)
        
        %% skip processing if cell was already visited or is below the minimum volume
        if (regionProps{i}(j).Area < minVolume || visitedIndices{i}(j))
            continue;
        end
        
        %% get the current object and add it to the result image
        currentFrame = i;
        currentObject = regionProps{i}(j).PixelIdxList;
        resultImages{currentFrame}(currentObject) = currentTrackingId;
        
        %% use transformed object, if registration correction is enabled
        if (useRegistrationBasedCorrection == true)
            currentObject = regionPropsTransformed{i}(j).PixelIdxList;
        end
        
        %% set visited status of the current object
        visitedIndices{i}(j) = 1;
        
        %% track current object as long as possible
        while currentFrame > 1
            
            %% identify the potential matches
            potentialMatches = segmentationImages{currentFrame-1}(currentObject);
            potentialMatchIndices = unique(potentialMatches);
            potentialMatchIndices(potentialMatchIndices == 0) = [];
            
            %% determine the best match
            matchCounts = zeros(size(potentialMatchIndices));
            diceIndices = zeros(size(potentialMatchIndices));
            for k=1:length(matchCounts)
                matchCounts(k) = sum(potentialMatches == potentialMatchIndices(k));
                diceIndices(k) = 2 * (matchCounts(k)) / (length(currentObject) + regionProps{currentFrame-1}(potentialMatchIndices(k)).Area);
            end
            
            %% set the match index as the maximum overlap segment
            [maxOverlap, maxIndex] = max(diceIndices);
            matchIndex = potentialMatchIndices(maxIndex);
            
            %% skip if the linked object is already taken or if no object was found
            if (isempty(matchIndex) || visitedIndices{currentFrame-1}(matchIndex) > 0)
                break;
            end
            
            %% update the matched object and add it to the next result image
            visitedIndices{currentFrame-1}(matchIndex) = 1;
            currentObject = regionProps{currentFrame-1}(matchIndex).PixelIdxList;
            resultImages{currentFrame-1}(currentObject) = currentTrackingId;
            
            %% use the transformed current object for better alignment
            if (useRegistrationBasedCorrection == true && length(regionPropsTransformed{currentFrame-1}) >= matchIndex)
                currentObject = regionPropsTransformed{currentFrame-1}(matchIndex).PixelIdxList;
            end
            
            %% decrement frame counter
            currentFrame = currentFrame - 1;
        end
        
        %% increase tracking id
        currentTrackingId = currentTrackingId + 1;
    end
    
    %% print status message
    fprintf('Finished tracking frame %i / %i\n', numFrames-i+1, numFrames);
end

%% save result images
parfor i=1:numFrames
    saveastiff(uint16(resultImages{i}), [outputFolderTracked strrep(inputFiles(i).name, '.tif', '_Tracked.tif')], options);
end