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

%% add dependencies
addpath('../ThirdParty/saveastiff_4.3/');

%% specify the input folders
tempFolder = tempdir;

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

%% setup input output directories
inputDir = uigetdir('X:\Data\KnautNYU_PrimordiumCellSegmentation\20210506_Stiching_EGFP_Caax-H2a-mCherry\Nuclei\', 'Select input folder containing the raw images of the nuclei.');
inputDir = [inputDir filesep];

%% get the output directory
outputRoot = uigetdir('C:\Users\stegmaier\Downloads\GlobalOutputTest\', 'Specify output directory for the identified transformations.');
outputRoot = [outputRoot filesep];

%% create transformations dir
transformationDir = [outputRoot 'Transformations' filesep];
if (~isfolder(transformationDir)); mkdir(transformationDir); end

registrationParameters = [pwd filesep 'parameters_BSpline.txt'];
if (~isfile(registrationParameters))
    uigetfile('*.txt', 'Select the transformation parameters file name.');
end

%% get the input files
inputFiles = dir([inputDir '*.tif']);
numFrames = length(inputFiles);

%% perform registration for all frames
for i=1:(numFrames-1)
    
    %% load the fixed (t-1) and the moving (t) image
    fixedImage = loadtiff([inputDir inputFiles(i).name]);
    movingImage = loadtiff([inputDir inputFiles(i+1).name]);
    
    %% normalize the input images and perform binary thresholding to primarily apply the transformation on the foreground
    normalizedFixedImage = double(fixedImage) / max(double(fixedImage(:)));
    normalizedMovingImage = double(movingImage) / max(double(movingImage(:)));
    maskThresholdFixed = graythresh(normalizedFixedImage);
    maskThresholdMoving = graythresh(normalizedMovingImage);
    fixedMaskImage = normalizedFixedImage >= maskThresholdFixed;
    movingMaskImage = normalizedMovingImage >= maskThresholdMoving;
    
    %% slightly enlarge the mask
    fixedMaskImage = imdilate(fixedMaskImage, strel('sphere', 2));
    movingMaskImage = imdilate(movingMaskImage, strel('sphere', 2));
    maskImage = max(fixedMaskImage, movingMaskImage);
    
    %% save temporary images for processing
    clear options;
    options.overwrite = true;
    saveastiff(uint16(fixedImage), [tempFolder 'fixedImage.tif'], options);
    saveastiff(uint16(movingImage), [tempFolder 'movingImage.tif'], options);
    saveastiff(uint16(maskImage), [tempFolder 'maskImageFixed.tif'], options);
    saveastiff(uint16(maskImage), [tempFolder 'maskImageMoving.tif'], options);
    
    %% perform the registration
    elastixCommand = [elastixPath 'elastix.sh ' ...
        '-f ' tempFolder 'fixedImage.tif ' ...
        '-m ' tempFolder 'movingImage.tif ' ...
        '-fMask ' tempFolder 'maskImageFixed.tif ' ...
        '-mMask ' tempFolder 'maskImageMoving.tif ' ...
        '-out ' tempFolder ' ' ...
        '-p ' registrationParameters];
    if (ispc)
        elastixCommand = strrep(elastixCommand, 'elastix.sh', 'elastix.exe');
    end
    system(elastixCommand);
    
    %% copy the result transformation
    tempResultFile = fopen([tempFolder 'TransformParameters.0.txt'], 'rb');
    resultFile = fopen(sprintf('%stransformation_t=%03d.txt', transformationDir, i+1), 'wb');
    
    %% replace interpolation and data type tags in the transformation file
    currentLine = fgetl(tempResultFile);
    while (~feof(tempResultFile))
        currentLine = strrep(currentLine, '(ResampleInterpolator "FinalBSplineInterpolator")', '(ResampleInterpolator "FinalNearestNeighborInterpolator")');
        currentLine = strrep(currentLine, '(ResultImageFormat "mhd")', '(ResultImageFormat "tif")');
        currentLine = strrep(currentLine, '(ResultImagePixelType "short")', '(ResultImagePixelType "unsigned short")');
        currentLine = strrep(currentLine, '(CompressResultImage "false")', '(CompressResultImage "true")');
        
        fprintf(resultFile, '%s\n', currentLine);
        currentLine = fgetl(tempResultFile);
    end
    fclose(tempResultFile);
    fclose(resultFile);
end