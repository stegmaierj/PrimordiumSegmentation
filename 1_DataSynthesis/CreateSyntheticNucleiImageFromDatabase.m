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

%% add library for fast tiff writing
addpath('../ThirdParty/saveastiff_4.3/');

%% load the nuclei data base
%load('nucleiDataBase.mat');
numTemplateImages = size(nucleiDataBase,1);
numTemplates = size(nucleiDataBase,2);

%%%%%%%% PARAMETERS %%%%%%%%
imageSize = [256, 256, 128];    %% image size of the created images
allowedOverlap = 10;            %% allowed overlap between neighboring cells
regionPadding = 4;              %% region padding used by the database creation script
safetyBorder = 20;              %% mask border that ensures no objects are positioned close to the boundaries
numImages = 210;                %% number of images that should be created
minNumObjects = 50;             %% minimum number of objects per image
additionalObjects = 200;        %% a random number <= this number is generated and added to the minimum number of objects

gaussianSigmaBefore = 1.0;            %% smoothing to be applied after the images were generated (used 0.4 for NYU data and 2.0 for CTC data)
gaussianWindowBefore = max(1, round(gaussianSigmaBefore)) * 3 + 1;  %% Gaussian window, defaulting to 3*sigma (i.e., covering 99% of the values)

gaussianSigmaAfter = 0.4;            %% smoothing to be applied after the images were generated (used 0.4 for NYU data and 1.0 for CTC data)
gaussianWindowAfter = max(1, round(gaussianSigmaAfter)) * 3 + 1;  %% Gaussian window, defaulting to 3*sigma (i.e., covering 99% of the values)

if (mod(gaussianWindowBefore, 2) == 0)
    gaussianWindowBefore = gaussianWindowBefore + 1;
end

if (mod(gaussianWindowAfter, 2) == 0)
    gaussianWindowAfter = gaussianWindowAfter + 1;
end

useGaussianBeforeNoise = true;
useGaussianAfterNoise = false;

clear settings;
settings.isotropicResultImage = false;
settings.zscale = 1;
settings.scaleFactor = 0.1;
settings.rotateX = false;
settings.rotateY = false;
settings.rotateZ = true;

settings.minCellDistance = 1;
settings.formClusters = true;
settings.clusterProbabilityThreshold = 0.9; %% random number between [0,1] is drawn. If value is smaller than the threshold, the object will be attached to a cluster, otherwise it will remain isolated.

settings.minZSlices = 3;
settings.meanStdImageDilation = 1;

settings.smoothMeanImage = true;
settings.meanImageSmoothingStdDev = 2;

%%%%%%%% PARAMETERS %%%%%%%%

%% get the output path
outputPath = uigetdir('X:\Projects\KnautNYU_PrimordiumCellSegmentation\2021_08_28_SyntheticNuclei\', 'Select output path for the generated images.');
outputPath = [outputPath filesep];
outputFolderLabelImages = [outputPath 'LabelImages' filesep];
outputFolderRawImages = [outputPath 'RawImages' filesep];
if (~isfolder(outputFolderLabelImages)); mkdir(outputFolderLabelImages); end
if (~isfolder(outputFolderRawImages)); mkdir(outputFolderRawImages); end

%% options for tiff writing
clear options;
options.compress = 'lzw';
options.overwrite = true;

%% was parfor before
parfor f=1:numImages
    
    %try
        %% per image settings
        numCells = minNumObjects + randperm(additionalObjects, 1);

        imageTemplate = randperm(numTemplateImages, 1);
        backgroundMean = nucleiDataBase(imageTemplate, 1).backgroundMean;
        backgroundStd = nucleiDataBase(imageTemplate, 1).backgroundStd;
        
        %% initialize label and statistics images
        labelImage = uint16(zeros(imageSize));      
        stdImage = backgroundStd * ones(imageSize);
        meanImage = backgroundMean * ones(imageSize);
        
        %% initialize the mask image and add the safety border to avoid sampling close to the border regions
        maskImage = uint16(zeros(imageSize));
        maskImage(:,:,1:safetyBorder) = 1;
        maskImage(:,1:safetyBorder,:) = 1;
        maskImage(1:safetyBorder,:,:) = 1;
        maskImage((end-safetyBorder):end,:,:) = 1;
        maskImage(:,(end-safetyBorder):end,:) = 1;
        maskImage(:,:,(end-safetyBorder):end) = 1;

        %% add the desired number of cells
        currentLabel = 1;
        for i=1:numCells

            %% pick random cell from data base
            currentTemplateId = randperm(numTemplates, 1);
            currentMaskSnippet = nucleiDataBase(imageTemplate, currentTemplateId).maskImage;
            currentMeanSnippet = nucleiDataBase(imageTemplate, currentTemplateId).meanImage;
            currentStdSnippet = nucleiDataBase(imageTemplate, currentTemplateId).stdImage;
            snippetSize = size(currentMaskSnippet);
            
            %% avoid degenerate cases
            if (snippetSize(3) < settings.minZSlices || sum(currentMaskSnippet(:)>0) == 0 || isempty(currentMaskSnippet))
                disp('Skipping Object... 1');
                continue;
            end
            
            if (settings.smoothMeanImage == true)
                currentMeanSnippet = imgaussfilt3(currentMeanSnippet, settings.meanImageSmoothingStdDev);
            end

            %% scale the images if anisotropy should be compensated
            if (settings.zscale ~= 1)
               currentMaskSnippet = imresize3(currentMaskSnippet, [snippetSize(1), snippetSize(2), round(snippetSize(3) * settings.zscale)], 'nearest');
               currentMeanSnippet = imresize3(currentMeanSnippet, [snippetSize(1), snippetSize(2), round(snippetSize(3) * settings.zscale)]);
               currentStdSnippet = imresize3(currentStdSnippet, [snippetSize(1), snippetSize(2), round(snippetSize(3) * settings.zscale)]);
            end

            %% prepare the mask and statistics snippets for augmentation
            currentMaskSnippet = imclose(currentMaskSnippet((regionPadding+1):(end-regionPadding), (regionPadding+1):(end-regionPadding), (regionPadding+1):(end-regionPadding)) > 0, strel('sphere', 2));
            currentMeanSnippet = currentMeanSnippet((regionPadding+1):(end-regionPadding), (regionPadding+1):(end-regionPadding), (regionPadding+1):(end-regionPadding)) .* double(imdilate(currentMaskSnippet > 0, strel('sphere', settings.meanStdImageDilation)));
            currentStdSnippet = currentStdSnippet((regionPadding+1):(end-regionPadding), (regionPadding+1):(end-regionPadding), (regionPadding+1):(end-regionPadding)) .* double(imdilate(currentMaskSnippet > 0, strel('sphere', settings.meanStdImageDilation)));

            %% avoid degenerate cases
            snippetSize = size(currentMaskSnippet);
            if (length(snippetSize) <= 2 || snippetSize(3) < settings.minZSlices || sum(currentMaskSnippet(:)>0) == 0 || isempty(currentMaskSnippet))
                disp('Skipping Object... 2');
                continue;
            end
            
            %% perform augmentation, i.e., randomly rotate the snippet in 3D and stretch it randomly along the different axes
            [currentMaskSnippet, currentMeanSnippet, currentStdSnippet] = PerformAugmentation(currentMaskSnippet, currentMeanSnippet, currentStdSnippet, settings);
            
            %% scale the images back to the anisotropic version should be compensated
            if (settings.zscale ~= 1 && settings.isotropicResultImage == false)
               snippetSize = size(currentMaskSnippet);
               currentMaskSnippet = imresize3(currentMaskSnippet, [snippetSize(1), snippetSize(2), round(snippetSize(3) / settings.zscale)], 'nearest');
               currentMeanSnippet = imresize3(currentMeanSnippet, [snippetSize(1), snippetSize(2), round(snippetSize(3) / settings.zscale)]);
               currentStdSnippet = imresize3(currentStdSnippet, [snippetSize(1), snippetSize(2), round(snippetSize(3) / settings.zscale)]);
            end

            %% avoid degenerate cases
            if (snippetSize(3) < settings.minZSlices || sum(currentMaskSnippet(:)>0) == 0)
                disp('Skipping Object... 2');
                continue;
            end
            
            %% compute the current snippet size and get the region properties
            snippetSize = size(currentMaskSnippet);
            regionProps = regionprops3(currentMaskSnippet, 'PrincipalAxisLength', 'Centroid', 'EquivDiameter', 'Volume');
            [~, maxRegionIndex] = max(regionProps.Volume);

            if (~isempty(maxRegionIndex))
                currentMaskSnippet = currentMaskSnippet == maxRegionIndex;
            else
                continue;
            end

            %% find a valid sampling location that is sufficiently distant to existing cells
            samplingLocation = FindValidSamplingLocation(maskImage, regionProps.EquivDiameter(maxRegionIndex) - allowedOverlap);
            
            %% find closest existing object
            objectIndices = find(labelImage(:) > 0);
            if (~isempty(objectIndices) && (settings.formClusters == true && settings.clusterProbabilityThreshold < rand))
                [xpos, ypos, zpos] = ind2sub(size(labelImage), objectIndices);
                
                currentDistances = sqrt((xpos - samplingLocation(1)).^2 + (ypos - samplingLocation(2)).^2 + (zpos - samplingLocation(3)).^2);
                [minDistance, minIndex] = min(currentDistances);
                
                closestObjectPixel = [xpos(minIndex), ypos(minIndex), zpos(minIndex)];
                
                [xposSnippet, yposSnippet, zposSnippet] = ind2sub(size(currentMaskSnippet), find(currentMaskSnippet(:) > 0));
                xposSnippet = xposSnippet + samplingLocation(1) - (size(currentMaskSnippet,1) / 2);
                yposSnippet = yposSnippet + samplingLocation(2) - (size(currentMaskSnippet,2) / 2);
                zposSnippet = zposSnippet + samplingLocation(3) - (size(currentMaskSnippet,3) / 2);
                
                currentSnippetDistances = sqrt((xposSnippet - closestObjectPixel(1)).^2 + (yposSnippet - closestObjectPixel(2)).^2 + (zposSnippet - closestObjectPixel(3)).^2);
                [minDistanceSnippet, minIndexSnippet] = min(currentSnippetDistances);
                closestSnippetPixel = [xposSnippet(minIndexSnippet), yposSnippet(minIndexSnippet), zposSnippet(minIndexSnippet)];
                
                movementDirection = (closestObjectPixel - closestSnippetPixel)';
                movementLength = norm(movementDirection);
                
                samplingLocation = samplingLocation + movementDirection * (movementLength + settings.minCellDistance);
                
                test = 1;
            end
            
            %% skip current cell if no valid sampling location was found
            if (isempty(samplingLocation))
                disp('No remaining position available for current min distance!');
                break;
            end

            %% determine region of interest based on the sampled point
            samplingLocation = round(samplingLocation - regionProps.Centroid(maxRegionIndex,:)');
            rangeX = samplingLocation(1):(samplingLocation(1) + snippetSize(1)-1);
            rangeY = samplingLocation(2):(samplingLocation(2) + snippetSize(2)-1);
            rangeZ = samplingLocation(3):(samplingLocation(3) + snippetSize(3)-1);

            if (rangeX(1) < 1 || rangeX(end) > size(labelImage,1) || ...
                    rangeY(1) < 1 || rangeY(end) > size(labelImage,2) || ...
                    rangeZ(1) < 1 || rangeZ(end) > size(labelImage,3))
                continue;
            end

            %% set the label image and the distribution images
            labelImage(rangeX, rangeY, rangeZ) = max(labelImage(rangeX, rangeY, rangeZ), uint16(currentMaskSnippet)*currentLabel);
%             
%             intensityWeight = (1 - 0.25 * (rand - 0.5));
%             
%             mean(currentMeanSnippet(currentMaskSnippet > 0))
%             
%             currentMeanSnippet = currentMeanSnippet - mean(currentMeanSnippet(currentMaskSnippet <= 0)) + backgroundMean;
            
            
            meanImage(rangeX, rangeY, rangeZ) = max(meanImage(rangeX, rangeY, rangeZ), intensityWeight * currentMeanSnippet); %% potentially use masked version with dilated nucleus mask.
            stdImage(rangeX, rangeY, rangeZ) = max(stdImage(rangeX, rangeY, rangeZ),  intensityWeight * currentStdSnippet); %% potentially use masked version with dilated nucleus mask.

            %% update the mask image used for sampling point computation
            maskImage = max(maskImage, uint16(labelImage > 0));
            currentLabel = currentLabel + 1;
        end

        %% sample current image and ensure that the limits fit
        if (useGaussianBeforeNoise == true)
            meanImage = imfilter(meanImage, fspecial3('gaussian', gaussianWindowBefore, gaussianSigmaBefore), 'replicate');
        end
        rawImage = meanImage + randn(imageSize) .* stdImage;
        if (useGaussianAfterNoise == true)
            rawImage = imfilter(rawImage, fspecial3('gaussian', gaussianWindowAfter, gaussianSigmaAfter), 'replicate');
        end
        rawImage(rawImage < 0) = 0;
        rawImage(rawImage > 65535) = 65535;

        %% save the current image
        saveastiff(uint16(labelImage), sprintf('%slabelImage%03d.tif', outputFolderLabelImages, f), options);
        saveastiff(uint16(rawImage), sprintf('%srawImage%03d.tif', outputFolderRawImages, f), options);

        %% status message
        disp(['Finished processing ' num2str(f) ' / ' num2str(numImages) ' ...']);
    
%     catch
%        disp(['Skipping ' num2str(f) ' due to error ...']); 
%     end
end