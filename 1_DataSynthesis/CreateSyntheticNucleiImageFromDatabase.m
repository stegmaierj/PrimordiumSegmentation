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
load('nucleiDataBase.mat');
numTemplates = length(nucleiDataBase);

%%%%%%%% PARAMETERS %%%%%%%%
imageSize = [256, 256, 128];    %% image size of the created images
allowedOverlap = 10;            %% allowed overlap between neighboring cells
regionPadding = 4;              %% region padding used by the database creation script
safetyBorder = 20;              %% mask border that ensures no objects are positioned close to the boundaries
numImages = 300;                %% number of images that should be created
minNumObjects = 50;             %% minimum number of objects per image
additionalObjects = 200;        %% a random number <= this number is generated and added to the minimum number of objects
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

parfor f=1:numImages
    
    %% per image settings
    numCells = minNumObjects + randperm(additionalObjects, 1);
    
    %% initialize label and statistics images
    labelImage = uint16(zeros(imageSize));
    stdImage = nucleiDataBase(randperm(numTemplates, 1)).backgroundStd * ones(imageSize);
    meanImage = nucleiDataBase(randperm(numTemplates, 1)).backgroundMean * ones(imageSize);
    
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
        currentMaskSnippet = nucleiDataBase(currentTemplateId).maskImage;
        currentMeanSnippet = nucleiDataBase(currentTemplateId).meanImage;
        currentStdSnippet = nucleiDataBase(currentTemplateId).stdImage;
        
        %% prepare the mask and statistics snippets for augmentation
        currentMaskSnippet = imclose(currentMaskSnippet((regionPadding+1):(end-regionPadding), (regionPadding+1):(end-regionPadding), (regionPadding+1):(end-regionPadding)) > 0, strel('sphere', 2));
        currentMeanSnippet = currentMeanSnippet((regionPadding+1):(end-regionPadding), (regionPadding+1):(end-regionPadding), (regionPadding+1):(end-regionPadding)) .* double(imdilate(currentMaskSnippet > 0, strel('sphere', 2)));
        currentStdSnippet = currentStdSnippet((regionPadding+1):(end-regionPadding), (regionPadding+1):(end-regionPadding), (regionPadding+1):(end-regionPadding)) .* double(imdilate(currentMaskSnippet > 0, strel('sphere', 2)));
        
        if (sum(currentMaskSnippet(:)>0) == 0)
            continue;
        end
        
        %% perform augmentation, i.e., randomly rotate the snippet in 3D and stretch it randomly along the different axes
        [currentMaskSnippet, currentMeanSnippet, currentStdSnippet] = PerformAugmentation(currentMaskSnippet, currentMeanSnippet, currentStdSnippet);
        
        %% compute the current snippet size and get the region properties
        snippetSize = size(currentMaskSnippet);
        regionProps = regionprops3(currentMaskSnippet, 'PrincipalAxisLength', 'Centroid', 'EquivDiameter', 'Volume');
        [~, maxRegionIndex] = max(regionProps.Volume);
        
        currentMaskSnippet = currentMaskSnippet == maxRegionIndex;
        
        %% find a valid sampling location that is sufficiently distant to existing cells
        samplingLocation = FindValidSamplingLocation(maskImage, regionProps.EquivDiameter(maxRegionIndex) - allowedOverlap);
        
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
        meanImage(rangeX, rangeY, rangeZ) = max(meanImage(rangeX, rangeY, rangeZ), (1 - 0.25 * (rand - 0.5)) * currentMeanSnippet); %% potentially use masked version with dilated nucleus mask.
        stdImage(rangeX, rangeY, rangeZ) = max(stdImage(rangeX, rangeY, rangeZ),  (1 - 0.25 * (rand - 0.5)) * currentStdSnippet); %% potentially use masked version with dilated nucleus mask.
        
        %% update the mask image used for sampling point computation
        maskImage = max(maskImage, uint16(labelImage > 0));
        currentLabel = currentLabel + 1;
    end
    
    %% sample current image and ensure that the limits fit
    rawImage = meanImage + randn(imageSize) .* stdImage;
    rawImage = imfilter(rawImage, fspecial3('gaussian', 3, 0.4));
    rawImage(rawImage < 0) = 0;
    rawImage(rawImage > 65535) = 65535;
    
    %% save the current image
    saveastiff(uint16(labelImage), sprintf('%slabelImage%03d.tif', outputFolderLabelImages, f), options);
    saveastiff(uint16(rawImage), sprintf('%srawImage%03d.tif', outputFolderRawImages, f), options);
    
    %% status message
    disp(['Finished processing ' num2str(f) ' / ' num2str(numImages) ' ...']);
end