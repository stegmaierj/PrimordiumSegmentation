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

%% get the input folders
nucleiFolder = uigetdir('/Users/jstegmaier/ScieboDrive/Projects/2021/KnautNYU_PrimordiumCellSegmentation/Data/SyntheticNuclei/Segmentation/', 'Please select the input folder containing 3D images with segmented nuclei.');
rawImageFolder = uigetdir('/Users/jstegmaier/ScieboDrive/Projects/2021/KnautNYU_PrimordiumCellSegmentation/Data/SyntheticNuclei/Raw/', 'Please select the raw image folder containing nuclei.');

%%%%%%%% PARAMETERS %%%%%%%%
%% binary threshold used for identifying the background.
%% potentially needs to be adjusted depending on the input image.
intensityThreshold = 800;

%% window size used for estimating mean and standard deviation
windowSize = 5;

%% padding in pixel used for cropping the nuclei templates
regionPadding = 4;
%%%%%%%% PARAMETERS %%%%%%%%

%% initialize the nuclei data base
nucleiDataBase = struct();
currentIndex = 1;

%% parse input directories for valid files
nucleiFiles = dir([nucleiFolder '*.tif']);
rawFiles = dir([rawImageFolder '*.tif']);

%% add all templates found in the segmentation images to the nucleus data base
for f = 1:length(nucleiFiles)

    %% load segmentation and raw images
    nucleiImage = loadtiff([nucleiFolder nucleiFiles(f).name]);
    rawImage = double(loadtiff([rawImageFolder rawFiles(f).name]));
    imageSize = size(rawImage);

    %% obtain the background image to assess the background characteristics
    maskImage = rawImage > intensityThreshold;

    %% compute std. and mean images
    stdImage = stdfilt(rawImage, true(windowSize, windowSize, windowSize));
    meanImage = imfilter(rawImage, fspecial3('average', [windowSize , windowSize, windowSize]));

    %% compute background statistics
    backgroundMean = mean(rawImage(~maskImage));
    backgroundStd = std(rawImage(~maskImage));
    
    %% extract the region props of the current image
    regionProps = regionprops(nucleiImage, 'Area', 'PixelIdxList', 'BoundingBox');
    
    %% add all valid nuclei to the database of potential objects
    for i=1:length(regionProps)

        %% skip empty entries
        if (regionProps(i).Area <= 0 || sum(nucleiImage(:) > 0) == 0)
            continue;
        end

        %% compute the distribution parameters for the current region
        currentMeans = meanImage(regionProps(i).PixelIdxList);
        currentStds = stdImage(regionProps(i).PixelIdxList);

        %% determine the AABB and crop the padded region
        aabb = regionProps(i).BoundingBox;
        rangeX = round(max(1, aabb(1)-regionPadding)):round(min(imageSize(2), aabb(1)+aabb(4)+regionPadding));
        rangeY = round(max(1, aabb(2)-regionPadding)):round(min(imageSize(1), aabb(2)+aabb(5)+regionPadding));
        rangeZ = round(max(1, aabb(3)-regionPadding)):round(min(imageSize(3), aabb(3)+aabb(6)+regionPadding));

        %% specify output images
        nucleiDataBase(currentIndex).rawImage = rawImage(rangeY, rangeX, rangeZ);
        nucleiDataBase(currentIndex).maskImage = imclose(uint8(nucleiImage(rangeY, rangeX, rangeZ) == i), strel('sphere', 2));
        nucleiDataBase(currentIndex).maskImage = imdilate(imgaussfilt3(nucleiDataBase(currentIndex).maskImage, 2) > 0.5, strel('sphere', 1));
        nucleiDataBase(currentIndex).stdImage = stdImage(rangeY, rangeX, rangeZ);
        nucleiDataBase(currentIndex).meanImage = meanImage(rangeY, rangeX, rangeZ);
        nucleiDataBase(currentIndex).backgroundMean = backgroundMean;
        nucleiDataBase(currentIndex).backgroundStd = backgroundStd;

        %% increment counter for next object and show status
        currentIndex = currentIndex + 1;
        disp(['Finished extracting ' num2str(currentIndex) ' / ' num2str(length(regionProps)) ' regions ...']);
    end
end

%% save the nucleus database
save('nucleiDataBase.mat', 'nucleiDataBase');