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
addpath('../ThirdParty/saveastiff_4.3');

%% specify input and output folders
inputFolder = uigetdir('X:\Data\KnautNYU_PrimordiumCellSegmentation\20210715_Stiching_EGFP_Caax-H2a-mCherry_crop\Membranes\', 'Specify the input folder containing the images to normalize in tiff format.');
inputFolder = [inputFolder filesep];
outputFolder = uigetdir('X:\Data\KnautNYU_PrimordiumCellSegmentation\20210715_Stiching_EGFP_Caax-H2a-mCherry_crop\Membranes_Normalized\', 'Specify the output path to write the normalized images to.');
outputFolder = [outputFolder filesep]

%% change the following value according to the quantile you want to use for normalization
quantileValue = 0.001;

%% get the input files
inputFiles = dir([inputFolder '*.tif']);
clear options;
options.compress = 'lzw';
options.overwrite = true;

%% process images in parallel and apply the normalization
parfor i=1:length(inputFiles)
    
    %% only process existing files
    if (exist([outputFolder inputFiles(i).name], 'file'))
        continue;
    end
    
    %% load current image
    currentImage = double(loadtiff([inputFolder inputFiles(i).name]));
    
    %% identify non-zero pixels and compute quantiles on the remianing pixels
    validIndices = currentImage(:) > 0;
    lowerQuantile = quantile(currentImage(validIndices), quantileValue);
    upperQuantile = quantile(currentImage(validIndices), 1.0-quantileValue);
    
    %% apply normalization to the current image
    currentImage = (currentImage - lowerQuantile) / (upperQuantile - lowerQuantile);
    currentImage(currentImage < 0) = 0;
    currentImage(currentImage > 1) = 1;
    
    %% save the result image.
    saveastiff(uint16(65535*currentImage), [outputFolder inputFiles(i).name], options);
end