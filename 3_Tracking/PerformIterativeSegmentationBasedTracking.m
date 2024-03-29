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
skipFrames = 10; %% the allowed number of frames to skip for tracking an object
maxCellSize = 4000; % 3000-4000 for Weiyis Project; 5.5599e+03 + 2*1.2410e+04; for Baptiste's project
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
inputDir = [outputRoot 'Nuclei_Segmented_Corrected' filesep];
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

outputFolderTrackedMaxProj = [outputRoot 'Nuclei_Tracked_MaxProj' filesep];
if (~isfolder(outputFolderTrackedMaxProj)); mkdir(outputFolderTrackedMaxProj); end

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
        segmentationImages{i} = uint16(loadtiff([inputDir inputFiles(i).name]));
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

maxLabel = length(regionProps{1});
d_orgs = zeros(maxLabel, numFrames, 5);

for i=1:numFrames
    for j=1:length(regionProps{i})
        if (regionProps{i}(j).Area > 0)
            d_orgs(j,i,:) = [j, regionProps{i}(j).Area, regionProps{i}(j).Centroid];
        end
    end    
end

%% plot statistics of the different frames

volumeStatistics = zeros(numFrames, 3);

for i=1:numFrames
    validIndices = find(d_orgs(:,i,1) > 0);

    volumeStatistics(i, 1) = mean(d_orgs(validIndices, i, 2));
    volumeStatistics(i, 2) = std(d_orgs(validIndices, i, 2));
    volumeStatistics(i, 3) = median(d_orgs(validIndices, i, 2));
    volumeStatistics(i, 4) = 1.4826*median(abs(median(d_orgs(validIndices, i, 2)) - (d_orgs(validIndices, i, 2))));

%     figure(2);
%     histogram(d_orgs(validIndices, i, 2));
%     title(['Frame: ' num2str(i) ', Mean: ' num2str(volumeStatistics(i,1)) ', Std.Dev.: ' num2str(volumeStatistics(i,2)) ', Median: ' num2str(volumeStatistics(i,3))]);
% 
%     test = 1;
end



%% clear previous tracking results
clear visitedIndices;
for i=1:numFrames
    visitedIndices{i} = zeros(size(regionProps{i},1), 1); %#ok<SAGROW>
    resultImages{i}(:) = 0;
end

%% initialize current tracking id
currentTrackingId = 1;

%% iterate over all frames in a backward fashion
for i=numFrames:-1:1

    %% process all cells contained in the current frame
    for j=1:size(regionProps{i},1)
        
        if (i==4 && (j == 178 || j == 179))
            test = 1;
        end

        %% skip processing if cell was already visited or is below the minimum volume
        if (regionProps{i}(j).Area < minVolume || visitedIndices{i}(j))
            continue;
        end

        if (currentTrackingId == 131)
            test = 1;
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
            
            matchIndices = [];
            matchDiceIndices = [];
            for s=1:skipFrames

                if ((currentFrame-s) < 1)
                    matchIndices = [];
                    break;
                end

                %% identify the potential matches
                currentIndex = unique(segmentationImages{currentFrame}(currentObject));
                potentialMatches = segmentationImages{currentFrame-s}(currentObject);
                potentialMatchIndices = unique(potentialMatches);
                potentialMatchIndices(potentialMatchIndices == 0) = [];

                if (isempty(potentialMatchIndices))
                    continue;
                end
                
                for m=potentialMatchIndices'

                    %% TODO: add forward consistency check for largest overlap matches!
                    matchObject = regionProps{currentFrame-s}(m).PixelIdxList;
                    potentialMatchesForward = segmentationImages{currentFrame}(matchObject);
                    potentialMatchIndicesFw = unique(potentialMatchesForward);
                    potentialMatchIndicesFw(potentialMatchIndicesFw == 0) = [];
        
                    matchCountsFw = zeros(size(potentialMatchIndicesFw));
                    diceIndicesFw = zeros(size(potentialMatchIndicesFw));
                    for k=1:length(matchCountsFw)
                        matchCountsFw(k) = sum(potentialMatchesForward == potentialMatchIndicesFw(k));
                        diceIndicesFw(k) = 2 * (matchCountsFw(k)) / (length(potentialMatchesForward) + regionProps{currentFrame}(potentialMatchIndicesFw(k)).Area);
                    end
        
                    [maxOverlapFw, maxIndexFw] = max(diceIndicesFw);
                    matchIndexFw = potentialMatchIndicesFw(maxIndexFw);
        
                    
                    %% skip if the linked object is already taken or if no object was found
                    if (ismember(matchIndexFw, currentIndex))
                        matchIndices = [matchIndices, m];
                        matchDiceIndices = [matchDiceIndices, maxOverlapFw];
                    end
                end
                

%                 %% determine the best match
%                 matchCounts = zeros(size(potentialMatchIndices));
%                 diceIndices = zeros(size(potentialMatchIndices));
%                 for k=1:length(matchCounts)
%                     matchCounts(k) = sum(potentialMatches == potentialMatchIndices(k));
%                     diceIndices(k) = 2 * (matchCounts(k)) / (length(currentObject) + regionProps{currentFrame-s}(potentialMatchIndices(k)).Area);
%                 end
% 
%                 [sortedDiceIndices, sortedIndices] = sort(diceIndices, 'descend');
% 
%                 if (length(sortedDiceIndices) > 1)
%                     secondNNRatio = sortedDiceIndices(2) / sortedDiceIndices(1);
%                     if (secondNNRatio > secondNNRatioThreshold)
%                         test = 1;
%                         disp('Oversegmentation Detected!');
%                     end
%                 else
% 
%                     %% set the match index as the maximum overlap segment
%                     [maxOverlap, maxIndex] = max(diceIndices);
%                     matchIndex = potentialMatchIndices(maxIndex);
%         
%                     %% skip if the linked object is already taken or if no object was found
%                     if (~isempty(matchIndex))
%                         currentSkipFrame = s;
%                         break;
%                     end
%                 end

                % skip if the linked object is already taken or if no object was found
                if (~isempty(potentialMatchIndices))
                    currentSkipFrame = s;
                    break;
                end
            end

            if (isempty(matchIndices))
                break;
            end

            %% add intermediate masks for frames where the segmentation is missing
            %% TODO: evenly distribute between current and previous frame.
            for s=1:(currentSkipFrame-1)
                resultImages{currentFrame-s}(currentObject) = currentTrackingId;
            end

            %% select only best match if the current segment is a potentially undersegmented object
            [sortedDiceIndices, sortedIndices] = sort(matchDiceIndices, 'descend');

            currentProbability = pdf('Normal', length(currentObject), volumeStatistics(currentFrame, 3), volumeStatistics(currentFrame, 4));

            if (length(currentObject) > maxCellSize)
                matchIndices = matchIndices(sortedIndices(1));
            end


            currentObject = [];
            for matchIndex=matchIndices

                if (visitedIndices{currentFrame-currentSkipFrame}(matchIndex) > 0)
                    continue;
                end

                %% update the matched object and add it to the next result image
                visitedIndices{currentFrame-currentSkipFrame}(matchIndex) = 1;
                currentObject = [currentObject; regionProps{currentFrame-currentSkipFrame}(matchIndex).PixelIdxList];
            end

            resultImages{currentFrame-currentSkipFrame}(currentObject) = currentTrackingId;
                
            %% TODO!!!! CHANGE !!! AFTER !!! ADDITION OF MULTIPLE SEGMENTS !!! use the transformed current object for better alignment
            if (useRegistrationBasedCorrection == true && length(regionPropsTransformed{currentFrame-currentSkipFrame}) >= matchIndex)
                currentObject = regionPropsTransformed{currentFrame-currentSkipFrame}(matchIndex).PixelIdxList;
            end
            
            %% decrement frame counter
            currentFrame = currentFrame - currentSkipFrame;
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
    saveastiff(uint16(max(resultImages{i}, [], 3)), [outputFolderTrackedMaxProj strrep(inputFiles(i).name, '.tif', '_TrackedMaxProj.tif')], options);
end