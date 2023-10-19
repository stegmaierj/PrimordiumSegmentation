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
exportCellSnippets = false;
preloadData = false; %% only set to false if debugging to avoid reloading the data at each call
minVolume = 100; %% excludes objects smaller than this volume
maxCellSize = 5.5599e+03 + 2*1.2410e+04; % 4000 for Weiyis Project;
%%%%%%%%%% PARAMETERS %%%%%%%%%%%

inputFolderRaw = 'X:/Data/KnautNYU_PrimordiumCellSegmentation/20211017_Stiching_EGFP_Caax-H2a-mCherry_crop/Nuclei/';
inputFilesRaw = dir([inputFolderRaw '*.tif']);

%% specify the output directory
outputRoot = uigetdir('C:\Users\stegmaier\Downloads\Processing\', 'Please select the output root directory ...');
outputRoot = [outputRoot filesep];
if (~isfolder(outputRoot)); mkdir(outputRoot); end

%% get the input directories
inputDir = [outputRoot 'Nuclei_Segmented' filesep];
if (~isfolder(inputDir))
    inputDir = uigetdir('I:\Projects\2021\KnautNYU_PrimordiumCellSegmentation\Processing\item_0012_GradientVectorFlowTrackingImageFilter\', 'Specify the input folder containing segmentation tiffs for each time point.');
    inputDir = [inputDir filesep];
else
   disp(['Using tracked images contained in the folder ' inputDir]);    
end

%% create the result folder for the tracking data
outputFolderTrackedCorrected = [outputRoot 'Nuclei_Segmented_Corrected' filesep];
if (~isfolder(outputFolderTrackedCorrected)); mkdir(outputFolderTrackedCorrected); end

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
    clear rawImages;
    clear resultImages;
    clear segmentationImages;
    clear regionProps;
    
    %% load and transform all input images
    parfor i=1:numFrames
        
        %% load the current image and obtain the region props
        rawImages{i} = uint16(loadtiff([inputFolderRaw inputFilesRaw(i).name]));
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

figure(3); clf; hold on;
plot(smooth(volumeStatistics(:,3), 0.1,'rloess'), '-r');
plot(smooth(volumeStatistics(:,3), 0.1,'rloess') + smooth(volumeStatistics(i, 4), 0.1,'rloess'), '--g');
plot(smooth(volumeStatistics(:,3), 0.1,'rloess') - smooth(volumeStatistics(i, 4), 0.1,'rloess'), '--g');

clear options;
options.overwrite = true;
options.compress = 'lzw';

for i=1:numFrames

    prevFrame = max(1, i-1);
    currFrame = i;
    nextFrame = min(i+1, numFrames);
    currentLabel = 1;

    for j=1:length(regionProps{i})

        if (regionProps{i}(j).Area <= 0)
            continue;
        end

        matchIndicesFw = [];
        matchDiceIndicesFw = [];
        matchIndicesBw = [];
        matchDiceIndicesBw = [];
        probFw = [];
        probBw = [];

        currentObject = regionProps{i}(j).PixelIdxList;
    
        %% identify the potential matches in the next frame
        currentIndex = unique(segmentationImages{currFrame}(currentObject));
        potentialMatches = segmentationImages{nextFrame}(currentObject);
        potentialMatchIndices = unique(potentialMatches);
        potentialMatchIndices(potentialMatchIndices == 0) = [];

        for m=potentialMatchIndices'

            if (isempty(m))
                continue;
            end

            %% TODO: add forward consistency check for largest overlap matches!
            matchObject = regionProps{nextFrame}(m).PixelIdxList;
            potentialMatchesForward = segmentationImages{currFrame}(matchObject);
            potentialMatchIndicesFw = unique(potentialMatchesForward);
            potentialMatchIndicesFw(potentialMatchIndicesFw == 0) = [];

            matchCountsFw = zeros(size(potentialMatchIndicesFw));
            diceIndicesFw = zeros(size(potentialMatchIndicesFw));
            for k=1:length(matchCountsFw)
                matchCountsFw(k) = sum(potentialMatchesForward == potentialMatchIndicesFw(k));
                diceIndicesFw(k) = 2 * (matchCountsFw(k)) / (length(potentialMatchesForward) + regionProps{currFrame}(potentialMatchIndicesFw(k)).Area);
            end

            [maxOverlapFw, maxIndexFw] = max(diceIndicesFw);
            matchIndexFw = potentialMatchIndicesFw(maxIndexFw);

            
            %% skip if the linked object is already taken or if no object was found
            if (ismember(matchIndexFw, currentIndex))
                matchIndicesFw = [matchIndicesFw, m];
                matchDiceIndicesFw = [matchDiceIndicesFw, maxOverlapFw];
                probFw = [probFw, pdf('Normal', length(matchObject), volumeStatistics(i,3), volumeStatistics(i,4))];
            end
        end

        %% identify the potential matches in the previous frame
        currentIndex = unique(segmentationImages{currFrame}(currentObject));
        potentialMatches = segmentationImages{prevFrame}(currentObject);
        potentialMatchIndices = unique(potentialMatches);
        potentialMatchIndices(potentialMatchIndices == 0) = [];

        for m=potentialMatchIndices'

            if (isempty(m))
                continue;
            end

            %% TODO: add forward consistency check for largest overlap matches!
            matchObject = regionProps{prevFrame}(m).PixelIdxList;
            potentialMatchesBackward = segmentationImages{currFrame}(matchObject);
            potentialMatchIndicesBw = unique(potentialMatchesBackward);
            potentialMatchIndicesBw(potentialMatchIndicesBw == 0) = [];

            matchCountsBw = zeros(size(potentialMatchIndicesBw));
            diceIndicesBw = zeros(size(potentialMatchIndicesBw));
            for k=1:length(matchCountsBw)
                matchCountsBw(k) = sum(potentialMatchesBackward == potentialMatchIndicesBw(k));
                diceIndicesBw(k) = 2 * (matchCountsBw(k)) / (length(potentialMatchesBackward) + regionProps{currFrame}(potentialMatchIndicesBw(k)).Area);
            end

            [maxOverlapBw, maxIndexBw] = max(diceIndicesBw);
            matchIndexBw = potentialMatchIndicesBw(maxIndexBw);

            
            %% skip if the linked object is already taken or if no object was found
            if (ismember(matchIndexBw, currentIndex))
                matchIndicesBw = [matchIndicesBw, m];
                matchDiceIndicesBw = [matchDiceIndicesBw, maxOverlapBw];
                probBw = [probBw, pdf('Normal', length(matchObject), volumeStatistics(i,3), volumeStatistics(i,4))];
            end
        end

        volumePrevObject = 0;
        for m=matchIndicesBw
            volumePrevObject = volumePrevObject + regionProps{prevFrame}(m).Area;
        end
        probPrevObject = pdf('Normal', volumePrevObject, volumeStatistics(i,3), volumeStatistics(i,4));

        volumeNextObject = 0;
        for m=matchIndicesFw
            volumeNextObject = volumeNextObject + regionProps{nextFrame}(m).Area;
        end
        probNextObject = pdf('Normal', volumeNextObject, volumeStatistics(i,3), volumeStatistics(i,4));
        probCurrObject = pdf('Normal', length(currentObject), volumeStatistics(i,3), volumeStatistics(i,4));

        if (length(matchIndicesBw) == 1 && length(matchIndicesFw) == 1)
        
            resultImages{i}(currentObject) = currentLabel;

            disp('Unambiguous match detected!');
%         elseif (length(matchIndicesBw) > 1 && max(probBw) > probCurrObject && ...
%                 length(matchIndicesFw) > 1 && max(probFw) > probCurrObject)
%             disp('Potential undersegmentation detected!');
%             %% TODO: handle over segmentation from both sides!
% 
% 
% 
%             resultImages{i}(currentObject) = 2;
%             test = 1;

        elseif (length(matchIndicesBw) > 1 && max(probBw) > probCurrObject)

            %% handle undersegmentation from previous frame
            disp('Potential undersegmentation detected!');

            [maxValue, maxIndex] = max(matchDiceIndicesBw);
            maxIndex = matchIndicesBw(maxIndex);

            potentialMatchesBackward = segmentationImages{prevFrame}(currentObject);
            potentialMatchIndicesBw = unique(potentialMatchesBackward);
            potentialMatchIndicesBw(potentialMatchIndicesBw == 0) = [];

            largestOverlap = currentObject(potentialMatchesBackward ~= maxIndex);

            resultImages{i}(currentObject) = currentLabel;
            
            currentLabel = currentLabel + 1;

            resultImages{i}(largestOverlap) = currentLabel;


        elseif (length(matchIndicesFw) > 1 && max(probFw) > probCurrObject)
            disp('Potential undersegmentation detected!');

            [maxValue, maxIndex] = max(matchDiceIndicesFw);
            maxIndex = matchIndicesFw(maxIndex);

            potentialMatchesForward = segmentationImages{nextFrame}(currentObject);
            potentialMatchIndicesFw = unique(potentialMatchesForward);
            potentialMatchIndicesFw(potentialMatchIndicesFw == 0) = [];

            largestOverlap = currentObject(potentialMatchesForward ~= maxIndex);

            resultImages{i}(currentObject) = currentLabel;
            
            currentLabel = currentLabel + 1;

            resultImages{i}(largestOverlap) = currentLabel;


        else
            %% other case
            resultImages{i}(currentObject) = currentLabel;

        end

        %% increment the current label
        currentLabel = currentLabel+1;
    end

    saveastiff(uint16(resultImages{i}), [outputFolderTrackedCorrected inputFiles(i).name], options);

    test = 1;
end

%%%% EXPORT SINGLE CELL SNIPPETS %%%%
if (exportCellSnippets == true)
    outputPath = 'C:/Users/stegmaier/Downloads/DebugVolume/';
    
    snippetWidth = 128;
    
    clear options;
    options.overwrite = true;
    options.compress = 'lzw';
    options.color = true;
    
    labelColor = [1,0,0];
    minIntensity = 0;
    
    for i=1:maxLabel
    
        validFrames = find(d_orgs(i,:,1) > 0);
    
        resultSnippet = uint16(zeros(snippetWidth,snippetWidth, 3, length(validFrames)));
        currentFrame = 1;
    
        mkdir(sprintf('%s/cell_%04d/', outputPath, i));
    
        for j=validFrames
    
            currentCentroid = regionProps{j}(i).Centroid;
            currentAABB = regionProps{j}(i).BoundingBox;
    
            width = round(currentAABB(4));
            height = round(currentAABB(5));
            depth = round(currentAABB(6));
    
            rangeX = max(1, round(currentCentroid(1) - snippetWidth/2):min(size(segmentationImages{j}, 2), round(currentCentroid(1) + snippetWidth/2)));
            rangeY = max(1, round(currentCentroid(2) - snippetWidth/2):min(size(segmentationImages{j}, 1), round(currentCentroid(2) + snippetWidth/2)));
            rangeZ = max(1, round(currentCentroid(3) - snippetWidth/2):min(size(segmentationImages{j}, 3), round(currentCentroid(3) + snippetWidth/2)));
    
            
            currentSnippet = uint16(segmentationImages{j}(rangeY, rangeX, rangeZ) == i);
            rawSnippet = rawImages{j}(rangeY, rangeX, rangeZ);
    
            resultSnippetR = imadjust(max(rawSnippet, [], 3));
            resultSnippetG = imadjust(max(rawSnippet, [], 3));
            resultSnippetB = imadjust(max(rawSnippet, [], 3));
    
            maxProjCurrentSnippet = max(currentSnippet, [], 3);
            currentIndices = find(maxProjCurrentSnippet(:) > 0);
    
            resultSnippetR(currentIndices) = max(minIntensity, resultSnippetR(currentIndices)) .* labelColor(1);
            resultSnippetG(currentIndices) = max(minIntensity, resultSnippetG(currentIndices)) .* labelColor(2);
            resultSnippetB(currentIndices) = max(minIntensity, resultSnippetB(currentIndices)) .* labelColor(3);
            
            resultSnippet = cat(3, resultSnippetR, resultSnippetG, resultSnippetG);
    
            imwrite(uint16(resultSnippet), sprintf('%s/cell_%04d/cell_%04d_t=%03d.png', outputPath, i, i, j));
            %saveastiff(uint16(resultSnippet), sprintf('%s/cell_%04d/cell_%04d_t=%03d.tif', outputPath, i, i, j), options);
    
            currentFrame = currentFrame + 1;
        end
    
    
        
    
        test = 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


colorMap = lines(maxLabel);

figure(2); clf;
for i=17

    validIndices = find(squeeze(d_orgs(i,:,1)) > 0);

    values1 = max( d_orgs(i,validIndices(1:(end-1)), 2), d_orgs(i,validIndices(2:end), 2));
    values2 = min( d_orgs(i,validIndices(1:(end-1)), 2), d_orgs(i,validIndices(2:end), 2));

    %plot3(d_orgs(i,validIndices, 3), d_orgs(i,validIndices, 4), d_orgs(i,validIndices, 2), '-k', 'Color', colorMap(i,:));
    plot(validIndices(2:end), values1 ./ values2, '-k', 'Color', colorMap(i,:)); hold on;

    medianVolume = median(d_orgs(i,validIndices,2));
    meanVolume = mean(d_orgs(i,validIndices,2));
    stdVolume = std(d_orgs(i,validIndices,2));
    lowerQuantile = quantile(d_orgs(i,validIndices,2), 0.0);
    upperQuantile = quantile(d_orgs(i,validIndices,2), 0.95);

    if (length(validIndices) <= 1)
        continue;
    end
% 
%     plot(validIndices, medianVolume * ones(1,length(validIndices)), '-k', 'LineWidth', 2);
%     plot(validIndices, upperQuantile * ones(1,length(validIndices)), '--k', 'LineWidth', 1);
%     plot(validIndices, lowerQuantile * ones(1,length(validIndices)), '--k', 'LineWidth', 1);
    
    title(['CellID: ' num2str(i) ', Median Volume: ' num2str(medianVolume) ', Std. Dev. Volume: ' num2str(stdVolume)]);
%     set(gca, 'YLim', [0, 6000]);
%     set(gca, 'XLim', [0, 338])

    hold off;
    %view(0,0);
    test = 1;
end
