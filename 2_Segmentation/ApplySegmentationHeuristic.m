%%
% CellCycleGAN.
% Copyright (C) 2021 D. Eschweiler, Weiyi Qian, H. Knaut, J. Stegmaier
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
% ...
%
%%

%% Performs segmentation based on registration and clustering.

%% load additional dependencies
addpath('ThirdParty/saveastiff_4.3');

%% specify the input and output directories
if (ismac)
    elastixPath = '/Users/jstegmaier/ScieboDrive/Projects/2021/KnautNYU_PrimordiumCellSegmentation/Software/ThirdParty/Elastix/macOS/bin/';
    inputDir = '/Users/jstegmaier/ScieboDrive/Projects/2021/KnautNYU_PrimordiumCellSegmentation/Processing/item_0017_TwangSegmentation/';
    transformationDir = '/Users/jstegmaier/ScieboDrive/Projects/2021/KnautNYU_PrimordiumCellSegmentation/Processing/Transformations/';
    resultDir = '/Users/jstegmaier/ScieboDrive/Projects/2021/KnautNYU_PrimordiumCellSegmentation/';
elseif (ispc)
    elastixPath = 'D:\ScieboDrive\Projects\2021\KnautNYU_PrimordiumCellSegmentation\Software\ThirdParty\Elastix\Windows\';
    %rawImageDir = 'I:\Projects\2021\KnautNYU_PrimordiumCellSegmentation\Data\';
    
    rawImageDir = 'X:\Data\KnautNYU_PrimordiumCellSegmentation\20210506_Stiching_EGFP_Caax-H2a-mCherry\Nuclei\';
    inputDir = 'I:\Projects\2021\KnautNYU_PrimordiumCellSegmentation\Processing\item_0012_GradientVectorFlowTrackingImageFilter\';
    transformationDir = 'D:\ScieboDrive\Projects\2021\KnautNYU_PrimordiumCellSegmentation\Processing\Transformations\';
    resultDir = 'C:\Users\stegmaier\Downloads\CorrectedResult\';
    if (~isfolder(resultDir)); mkdir(resultDir); end
else
    disp('Linux not supported yet ...');
end
useRegistrationBasedCorrection = false;



%% get information about the input/output files
rawImageFiles = dir([rawImageDir '*.tif']);
inputFiles = dir([inputDir '*.tif']);
transformationFiles = dir([transformationDir '*.txt']);
numFrames = length(inputFiles);

imageSize = [2474, 176, 77];

minVolume = 500;

currentTrackId = 1;

%% perform tracking for all frames starting at the last one
for i=2:(numFrames-1)
    
    %% compute the region props
    currRawImage = loadtiff([rawImageDir rawImageFiles(i).name]);
    prevSegmentation = loadtiff([inputDir inputFiles(i-1).name]);
    currSegmentation = loadtiff([inputDir inputFiles(i).name]);
    nextSegmentation = loadtiff([inputDir inputFiles(i+1).name]);
    
    currSegmentationCorrected = currSegmentation;
    
    %     prevRawImage = loadtiff([rawImageDir rawImageFiles(i-1).name]);
    %     currRawImage = loadtiff([rawImageDir rawImageFiles(i).name]);
    %     nextRawImage = loadtiff([rawImageDir rawImageFiles(i+1).name]);
    
    prevRegionProps = regionprops(prevSegmentation, 'Area', 'Centroid', 'PixelIdxList', 'BoundingBox');
    currRegionProps = regionprops(currSegmentation, 'Area', 'Centroid', 'PixelIdxList', 'BoundingBox');
    nextRegionProps = regionprops(nextSegmentation, 'Area', 'Centroid', 'PixelIdxList', 'BoundingBox');
    
    nextNewLabel = max(currSegmentation(:)) + 1;
    
    for j=1:length(currRegionProps)
        
        currentId = j;
        currentVolume = currRegionProps(currentId).Area;
        
        prevMatches = prevSegmentation(currRegionProps(currentId).PixelIdxList);
        nextMatches = nextSegmentation(currRegionProps(currentId).PixelIdxList);
        
        prevLabels = unique(prevMatches);
        nextLabels = unique(nextMatches);
        prevLabels(prevLabels == 0) = [];
        nextLabels(nextLabels == 0) = [];
        
        prev2CurrMatches = zeros(length(prevLabels), 3);
        next2CurrMatches = zeros(length(nextLabels), 3);
        
        %% identify best match in current frame
        for k=1:length(prevLabels)
            [prev2CurrMatches(k,1), prev2CurrMatches(k,2)] = FindMaxIndex(currSegmentation(prevRegionProps(prevLabels(k)).PixelIdxList));
            prev2CurrMatches(k,3) = 2 * (prev2CurrMatches(k,2)) / (currentVolume + prevRegionProps(prevLabels(k)).Area);
        end
        
        for k=1:length(nextLabels)
            [next2CurrMatches(k,1), next2CurrMatches(k,2)] = FindMaxIndex(currSegmentation(nextRegionProps(nextLabels(k)).PixelIdxList));
            next2CurrMatches(k,3) = 2 * (next2CurrMatches(k,2)) / (currentVolume + nextRegionProps(nextLabels(k)).Area);
        end
        
        validMatchesPrev2Curr = find(prev2CurrMatches(:,1) == currentId);
        validMatchesNext2Curr = find(next2CurrMatches(:,1) == currentId);
        
        if (length(validMatchesPrev2Curr) > 1 && length(validMatchesNext2Curr) > 1)
            
            
            
            unionPrev = 0;
            for k=validMatchesPrev2Curr'
                unionPrev = unionPrev + prevRegionProps(prevLabels(k)).Area;
            end
            
            prevDiceBefore= sort(prev2CurrMatches(validMatchesPrev2Curr,3), 'descend');
            nextDiceBefore= sort(next2CurrMatches(validMatchesNext2Curr,3), 'descend');
            if ((prevDiceBefore(2) / prevDiceBefore(1)) > 0.5 && (nextDiceBefore(2) / nextDiceBefore(1)) > 0.5)
                disp('largest overlaps are almost identical');
                
                
                
                
                diceAfter = 2*(sum(prev2CurrMatches(validMatchesPrev2Curr,2))) / (currentVolume + unionPrev)
                disp('Potential under segmenation spotted!');
                
                if (diceAfter < 0.5)
                    continue;
                end               
                
                
                figure(1); clf;
                subplot(2,3,1);
                currAABB = round(currRegionProps(j).BoundingBox);
                rangeY = max(1, currAABB(1)):min(size(currSegmentation,2), (currAABB(1)+currAABB(4)-1));
                rangeX = max(1, currAABB(2)):min(size(currSegmentation,1), (currAABB(2)+currAABB(5)-1));
                rangeZ = max(1, currAABB(3)):min(size(currSegmentation,3), (currAABB(3)+currAABB(6))); %round(currRegionProps(j).Centroid(3)); %
                imagesc(max(prevSegmentation(rangeX, rangeY, rangeZ), [], 3));
                
                subplot(2,3,2);
                imagesc(max(currSegmentation(rangeX, rangeY, rangeZ), [], 3)); hold on;
                plot(currAABB(4)/2, currAABB(5)/2, '*r');
                
                subplot(2,3,3);
                imagesc(max(nextSegmentation(rangeX, rangeY, rangeZ), [], 3));
                title(['Dice After: ' num2str(diceAfter)]);
                
                
                subplot(2,3,5);
                imagesc(max(currRawImage(rangeX, rangeY, rangeZ), [], 3));
                
                subplot(2,3,6);
                
                %% perform the split
                currentRegion = double(currRawImage(rangeX, rangeY, rangeZ)) .* double(currSegmentation(rangeX, rangeY, rangeZ) == currentId);
                currentRegion = max(currentRegion(:)) - currentRegion;
                
                seedCandidates = nextLabels(validMatchesNext2Curr);
                
                labelImage = zeros(size(currentRegion));
                for k=seedCandidates'
                    
                    seedImage = (double(nextSegmentation(rangeX, rangeY, rangeZ)) .* (double(currSegmentation(rangeX, rangeY, rangeZ)) == currentId)) == k;
                    [xpos, ypos, zpos] = ind2sub(size(seedImage), find(seedImage(:) > 0));
                    
                    labelImage(round(mean(xpos)), round(mean(ypos)), round(mean(zpos))) = 1;
                end
                currentRegion = imimposemin(currentRegion, labelImage);
                
                currenRegionWs = watershed(currentRegion);
                splitImage = double(currenRegionWs > 0) .* double(currSegmentation(rangeX, rangeY, rangeZ) == currentId);
                imagesc(splitImage(:,:,round(size(splitImage,3)/2)));
                
                %% replace the original segment with the split one
                %% TODO
                oldRegionContent = currSegmentation(rangeX, rangeY, rangeZ);
                oldRegionContent(oldRegionContent == currentId) = 0;
                
                localRegionProps = regionprops(splitImage, 'PixelIdxList');
                for k=1:length(localRegionProps)
                    if (k==1)
                        oldRegionContent(localRegionProps(k).PixelIdxList) = currentId;
                    else
                        oldRegionContent(localRegionProps(k).PixelIdxList) = nextNewLabel;
                        nextNewLabel = nextNewLabel+1;
                    end
                end
                
                currSegmentationCorrected(rangeX, rangeY, rangeZ) = oldRegionContent;
            end
        end
        

        
        %% correct under segmentation errors
        
        %% correct over segmentation errors
        
        
    end
    
    clear options;
    options.overwrite = true;
    options.compress = 'lzw';

    saveastiff(currSegmentationCorrected, [resultDir inputFiles(i).name], options);


    test = 1;
    
    fprintf('Finished tracking of %i / %i frames ...\n', numFrames-i+1, numFrames);
end