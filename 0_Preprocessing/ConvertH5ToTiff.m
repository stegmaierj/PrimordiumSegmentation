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

%% add lib for fast tiff writing
addpath('../ThirdParty/saveastiff_4.3/');

%% specify the input file
[inputFile, inputFolder] = uigetfile('*.h5', 'Select input file in H5 format'); %% Example: D:\20210715_Stiching_EGFP Caax-H2a-mCherry_crop.h5
inputFile = [inputFolder inputFile];
inputInfo = h5info(inputFile);
[path, fileName, ext] = fileparts(inputFile);

%% specify the output path
outputRoot = uigetdir('C:\Users\stegmaier\Downloads\GlobalOutputTest\', 'Please specify the output directory for the extracted frames');
outputRoot = [outputRoot filesep];
outputChannel1 = [outputRoot 'Membranes\'];
outputChannel2 = [outputRoot 'Nuclei\'];
if (~isfolder(outputChannel1)); mkdir(outputChannel1); end
if (~isfolder(outputChannel2)); mkdir(outputChannel2); end

%% get the number of time points and channels
numTimePoints = length(inputInfo.Groups)-1;
numChannels = 1;

%% options for tiff writing
clear options;
options.overwrite = true;
options.compress = 'lzw';

%% iteratively extract all time points from the h5 file
for t=1:numTimePoints
    for c=1:numChannels
        
        %% get the current frame
        currentImage = h5read(inputFile, sprintf('/t%03d/channel%d', t-1, c-1));
        currentImage = permute(currentImage, [2, 1, 3]);
        
        %% handle channels differently
        if (c==1)
            outputPath = sprintf('%s%s_t=%03d_c=%02d.tif', outputChannel1, fileName, t-1, c-1);
        else
            outputPath = sprintf('%s%s_t=%03d_c=%02d.tif', outputChannel2, fileName, t-1, c-1);
        end
        
        %% save result image of current frame
        saveastiff(currentImage, outputPath, options);
    end
end