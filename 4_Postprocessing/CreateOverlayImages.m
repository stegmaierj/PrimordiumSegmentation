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

%% add additional headers
addpath('../ThirdParty/saveastiff_4.3/');

%% select input folders
rawInputFolder = uigetdir('', 'Please select the folder containing the raw images.');
segInputFolder = uigetdir('', 'Please select the folder containing the segmented images.');
if (rawInputFolder(end) ~= filesep); rawInputFolder = [rawInputFolder filesep]; end
if (segInputFolder(end) ~= filesep); segInputFolder = [segInputFolder filesep]; end
segInputFiles = dir([segInputFolder '*.tif']);
rawInputFiles = dir([rawInputFolder '*.tif']);

%% select output folder
outputFolder = uigetdir('', 'Please select the output folder.');
if (outputFolder(end) ~= filesep); outputFolder = [outputFolder filesep]; end
if (~isfolder(outputFolder)); mkdir(outputFolder); end

%% perform validity check for input folders
if ((length(segInputFolder) == 1 && segInputFolder == 0) || ...
    (length(rawInputFolder) == 1 && rawInputFolder == 0) || ...
    isempty(segInputFiles) || isempty(rawInputFiles) || ...
    length(rawInputFiles) ~= length(segInputFiles))
    disp('Invalid file selection. Please select valid input paths for the raw images and segmentation images and make sure the number of images matches.');
    return;
end

%% load the first frame assuming due to the backwards tracking the largest label is contained in this frame.
firstFrameImage = loadtiff([segInputFolder segInputFiles(1).name]);
maxLabel = max(firstFrameImage(:));
labelColors = lines(maxLabel);
labelColors = labelColors(randperm(maxLabel)', :);
labelColors = rand(size(labelColors));

%% specify rendering parameters
minIntensity = 0.5;
maxArea = 1000;

%% process all images in parallel
parfor f=1:length(segInputFiles)
        
    %% load the current images
    segmentationImage = loadtiff([segInputFolder segInputFiles(f).name]);
    rawImage = loadtiff([rawInputFolder rawInputFiles(f).name]);
    segmentationImage = max(segmentationImage, [], 3);
    rawImage = double(max(rawImage, [], 3));
    
    %% perform intensity normalization
    rawImage = (rawImage - quantile(rawImage(:), 0.001)) / (quantile(rawImage(:), 0.999) - quantile(rawImage(:), 0.001));
    rawImage(rawImage > 1) = 1;
    rawImage(rawImage < 0) = 0;
    
    %% identify the current region props    
    regionProps = regionprops(segmentationImage, 'Area', 'PixelIdxList');
    
    %% initialize the result image channels
    resultImageR = rawImage;
    resultImageG = rawImage;
    resultImageB = rawImage;
    
    %% color all regions according to their track color
    for i = 1:length(regionProps)
        
        %% skip large objects (used for background suppression)
        if (regionProps(i).Area > maxArea)
            continue;
        end
        
        %% set the color of the current object
        resultImageR(regionProps(i).PixelIdxList) = max(minIntensity, rawImage(regionProps(i).PixelIdxList)) * labelColors(i,1);
        resultImageG(regionProps(i).PixelIdxList) = max(minIntensity, rawImage(regionProps(i).PixelIdxList)) * labelColors(i,2);
        resultImageB(regionProps(i).PixelIdxList) = max(minIntensity, rawImage(regionProps(i).PixelIdxList)) * labelColors(i,3);
    end
    
    %% concatenate the color channels to the final result image
    resultImage = cat(3, resultImageR, resultImageG, resultImageB);
    
    %% write result image
    imwrite(resultImage, [outputFolder strrep(rawInputFiles(f).name, '.tif', '_MaxProj.png')]);
end