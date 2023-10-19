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

%%%%%%%%%%%%% USER-DEFINED PARAMETERS %%%%%%%%%%%%
useRegistrationBasedCorrection = true;          %% if enabled, segmentation results are backpropagated via elastic registration results
zscale = 2.5;                                   %% ratio of z vs. xy spacing to compute distances correctly
clusterCutoff = 5;                              %% objects closer than this radius are merged
numGradientSteps = 20;                           %% number of gradient descent steps made to center detections in their basins
%%%%%%%%%%%%% USER-DEFINED PARAMETERS %%%%%%%%%%%%

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
outputRoot = uigetdir('C:\Users\stegmaier\Downloads\ManualInitTest\', 'Please select the output root directory ...');
outputRoot = [outputRoot filesep];
if (~isfolder(outputRoot)); mkdir(outputRoot); end

%% create the result folder for the tracking data
outputFolderTracked = [outputRoot 'Nuclei_Tracked' filesep];
if (~isfolder(outputFolderTracked)); mkdir(outputFolderTracked); end

%% create the result folder for the tracking data maximum projections
outputFolderTrackedMaxProj = [outputRoot 'Nuclei_Tracked_MaxProj' filesep];
if (~isfolder(outputFolderTrackedMaxProj)); mkdir(outputFolderTrackedMaxProj); end

%% get the input directories
inputDirNucleiBasins = [outputRoot 'Nuclei_BasinsRaw' filesep];
if (~isfolder(inputDirNucleiBasins))
    inputDirNucleiBasins = uigetdir('I:\Projects\2021\KnautNYU_PrimordiumCellSegmentation\Processing\item_0012_GradientVectorFlowTrackingImageFilter\', 'Specify the input folder containing raw nuclei basins for each time point.');
    inputDirNucleiBasins = [inputDirNucleiBasins filesep];
else
    disp(['Using segmented images contained in the folder ' inputDirNucleiBasins]);
end

%% get the input directories
inputDirNucleiBasins3Class = [outputRoot 'Nuclei_Basins3Class' filesep];
if (~isfolder(inputDirNucleiBasins3Class))
    inputDirNucleiBasins3Class = uigetdir('I:\Projects\2021\KnautNYU_PrimordiumCellSegmentation\Processing\item_0012_GradientVectorFlowTrackingImageFilter\', 'Specify the input folder containing 3 class nuclei basins for each time point.');
    inputDirNucleiBasins3Class = [inputDirNucleiBasins3Class filesep];
else
    disp(['Using segmented images contained in the folder ' inputDirNucleiBasins3Class]);
end

%% get information about the input/output files
inputFilesNucleiBasins = dir([inputDirNucleiBasins '*.tif']);
inputFilesNucleiBasins3Class = dir([inputDirNucleiBasins3Class '*.tif']);
numFrames = length(inputFilesNucleiBasins);

%% get the input directories
inputDirNucleiSeg = [outputRoot 'Nuclei_Segmented' filesep];
if (~isfolder(inputDirNucleiBasins))
    inputDirNucleiSeg = uigetdir('I:\Projects\2021\KnautNYU_PrimordiumCellSegmentation\Processing\item_0012_GradientVectorFlowTrackingImageFilter\', 'Specify the input folder containing segmentation tiffs for each time point.');
    inputDirNucleiSeg = [inputDirNucleiSeg filesep];
else
    disp(['Using segmented images contained in the folder ' inputDirNucleiSeg]);
end

%% get information about the input/output files
inputFilesNucleiSeg = dir([inputDirNucleiSeg '*.tif']);

trackInitFile = [outputRoot 'Tracking_Initialization.tif'];
if (~isfile(trackInitFile))
    [trackInitFile, trackInitFilePath] = uigetfile('*.tif', 'Specify the segmentation image of the last frame to use for tracking initialization.');
    trackInitFile = [trackInitFilePath trackInitFile];
else
    disp(['Using the following image to initialize the tracking: ' trackInitFile]);
end

%% load the track initialization file
trackingInitialization = loadtiff(trackInitFile);

%% specify save options
clear options;
options.overwrite = true;
options.compress = 'lzw';

%% get the transformations directory
if (useRegistrationBasedCorrection == true)
    transformationDir = [outputRoot 'Transformations' filesep];
    if (~isfolder(transformationDir))
        transformationDir = uigetdir('D:\ScieboDrive\Projects\2021\KnautNYU_PrimordiumCellSegmentation\Processing\Transformations\', 'Specify the transformation directory');
        transformationDir = [transformationDir filesep];
    else
        disp(['Using transformations contained in the folder ' transformationDir]);
    end

    transformationFiles = dir([transformationDir '*.txt']);
end

%% initialize d_orgs
d_orgs = zeros(0, numFrames, 6);
trackingIdIndex = 6;

%% perform tracking backwards in time
for t=numFrames:-1:1

    %% load the current basin image
    currImage3Class = loadtiff([inputDirNucleiBasins3Class inputFilesNucleiBasins3Class(t).name]);
    currImage3Class = imgaussfilt3(currImage3Class, [2.5, 2.5, 1]);

    currImage = loadtiff([inputDirNucleiBasins inputFilesNucleiBasins(t).name]);
    currImage = currImage / quantile(abs(currImage(:)), 0.99);
    currImage = currImage + currImage3Class;

    [currGradX, currGradY, currGradZ] = imgradientxyz(currImage);

    %% initialize the result image
    resultImages{t} = zeros(size(currImage));

    %% load segmentation initialization
    if (t == numFrames && ~isempty(trackingInitialization))
        segImage = trackingInitialization;
    else
        segImage = loadtiff([inputDirNucleiSeg inputFilesNucleiSeg(t).name]);

        %% remove objects left right of the largest connected component
        segImageMaxProj = max(segImage, [], 3);
        segImageMaxProjBW = imgaussfilt(segImageMaxProj, 5) > 0;

        %% find largest connected component
        regionProps = regionprops(segImageMaxProjBW, 'Area', 'BoundingBox');
        maxIndex = 1;
        maxArea = 0;
        for j=1:length(regionProps)
            if (regionProps(j).Area > maxArea)
                maxArea = regionProps(j).Area;
                maxIndex = j;
            end
        end

        %% remove all detections that are farther right of the largest connected component
        segImage(:, round(regionProps(maxIndex).BoundingBox(1)+regionProps(maxIndex).BoundingBox(3)):end, :) = 0;
    end

    %% compute region props of the current image
    regionProps = regionprops(segImage, 'Area', 'Centroid');

    %% identify the maximum index to know where to continue with the track labels
    currentMaxIndex = 0;
    if (~isempty(d_orgs))
        currentMaxIndex = max(find(d_orgs(:,t,1)));
    end

    %% add the detections of the current frame to d_orgs
    currentLabel = 1;
    for j=1:length(regionProps)
        if (regionProps(j).Area > 0)
            d_orgs(currentMaxIndex+currentLabel,t,:) = [currentMaxIndex+currentLabel, regionProps(j).Area, round(regionProps(j).Centroid([2,1,3])), 0];
            currentLabel = currentLabel + 1;
        end
    end

    %% perform gradient descent based movement of all detections
    validIndices = find(d_orgs(:,t,1) > 0);

    for i=validIndices'

        %% get the current position and check for validity
        currPosition = squeeze(d_orgs(i, t, 3:5));

        if (sum(currPosition == 0))
            continue;
        end

        if (t < numFrames)
            %% perform the desired number of gradient move steps
            for j=1:numGradientSteps
    
                %% skip if position is invalid
                if (sum(currPosition == 0))
                    break;
                end
    
                %% get the gradient for the current position
                currGradient = [currGradY(currPosition(1), currPosition(2), currPosition(3)), ...
                    currGradX(currPosition(1), currPosition(2), currPosition(3)), ...
                    currGradZ(currPosition(1), currPosition(2), currPosition(3))]';
    
                %% perform gradient step using a discretized gradient
                currPosition = round(currPosition - sign(currGradient));
    
                %% ensure the position remains within the image boundaries
                currPosition(1) = max(1, min(currPosition(1), size(currImage, 1)));
                currPosition(2) = max(1, min(currPosition(2), size(currImage, 2)));
                currPosition(3) = max(1, min(currPosition(3), size(currImage, 3)));
            end
        end

        %% skip processing if position is invalid
        if (sum(currPosition == 0))
            continue;
        end

        %% save the updated position
        d_orgs(i, t, 3:5) = currPosition;

    end

    %% get all candidate positions for the current frame
    currentPositions = squeeze(d_orgs(validIndices, t, 3:5));
    currentPositions(:,3) = currentPositions(:,3) * zscale;

    %% perform clustering
    Z = linkage(currentPositions, 'ward');
    c = cluster(Z, 'Cutoff', clusterCutoff, 'Criterion','distance');

    %% TODO: add valid cells to the result image
    clusterIndices = unique(c);
    currentLabel = max(find(d_orgs(:, t, trackingIdIndex) > 0)) + 1;
    if (isempty(currentLabel))
        currentLabel = 1;
    end

    %% initialize temporary detections for the current frame
    d_orgs_temp = zeros(size(squeeze(d_orgs(:,t,:))));

    %% process all identified clusters
    for i=clusterIndices'

        %% get indices of all objects contained in the cluster and the ones that already have a track id
        currentObjects = validIndices(c == i);
        currentTrackedObjects = validIndices(c == i & d_orgs(validIndices,t,trackingIdIndex) > 0);

        %% combine positions of all cluster objects by averaging
        if (length(currentObjects) > 1)
            meanPosition = round(mean(squeeze(d_orgs(currentObjects, t, 3:5)), 1));
        else
            meanPosition = squeeze(d_orgs(currentObjects, t, 3:5));
        end

        %% add the object with a new tracking id or keep an existing one if available
        if (isempty(currentTrackedObjects))
            currentTrackingId = currentLabel;
            currentLabel = currentLabel + 1;
        else
            currentTrackingId = currentTrackedObjects(1);
        end

        %% update the temporary d_orgs variable
        d_orgs_temp(currentTrackingId, 1) = currentTrackingId;
        d_orgs_temp(currentTrackingId, 2) = 1;
        d_orgs_temp(currentTrackingId, 3:5) = meanPosition;
        d_orgs_temp(currentTrackingId, trackingIdIndex) = 1;

        %% add the current object to the result image
        resultImages{t}(meanPosition(1), meanPosition(2), meanPosition(3)) = currentTrackingId;
    end

    %% copy temporary results to d_orgs
    d_orgs(:, t, :) = d_orgs_temp;

    %% dilate the current objects to avoid information loss due to the elastic deformation
    resultImages{t} = imdilate(resultImages{t}, strel('sphere', 3));

    %% transform current result image to the previous frame
    if (useRegistrationBasedCorrection == true)

        %% save temporary image
        outputDirTemp = [tempdir num2str(i) filesep];
        inputFile = [outputDirTemp 'currImage.tif'];
        saveastiff(uint16(resultImages{t}), inputFile, options);

        %% apply transformation to the segmentation image
        transformixCommand = [elastixPath 'transformix.sh ' ...
            '-in ' inputFile ' ' ...
            '-out ' outputDirTemp ' ' ...
            '-tp ' transformationDir transformationFiles(t).name];

        if (ispc)
            transformixCommand = strrep(transformixCommand, 'transformix.sh', 'transformix.exe');
        end
        system(transformixCommand);

        %% load the transformed image at time point t and the previous image at t-1
        transformedCurrentImage = loadtiff([outputDirTemp 'result.tif']);

    else
        transformedCurrentImage = resultImages{t};
    end

    %% save the regionprops for later computations
    regionPropsTransformed = regionprops(transformedCurrentImage, 'Area', 'Centroid', 'PixelIdxList', 'BoundingBox');

    %% set the transformed positions as the intialization for the previous frame
    if (t > 1)
        for i=1:size(regionPropsTransformed)
            if (regionPropsTransformed(i).Area <= 0)
                d_orgs(i, t-1, :) = 0;
            else
                d_orgs(i, t-1, :) = d_orgs(i, t, :);
                d_orgs(i, t-1, trackingIdIndex) = 1;
                d_orgs(i, t-1, 3:5) = round(regionPropsTransformed(i).Centroid([2,1,3]));
            end
        end
    end

    %% save the current result images (note that a filtered version will later overwrite these results)
    saveastiff(uint16(resultImages{t}), [outputFolderTracked strrep(inputFilesNucleiBasins(t).name, '.tif', '_Tracked.tif')], options);
    saveastiff(uint16(max(resultImages{t}, [], 3)), [outputFolderTrackedMaxProj strrep(inputFilesNucleiBasins(t).name, '.tif', '_TrackedMaxProj.tif')], options);

    %% print status message
    disp(['Finished tracking ' num2str(numFrames - t) ' / ' num2str(numFrames)]);
end

%% analyze track lengths
trackLengths = zeros(size(d_orgs,1), 1);
endTimePoints = zeros(size(d_orgs,1), 1);
trackLengthThreshold = 25;
trackLabels = zeros(size(d_orgs,1), 1);
currentTrackLabel = 1;
for i=1:size(d_orgs,1)
    currentIndices = find(d_orgs(i, :, 1) > 0);
    trackLengths(i) = length(currentIndices);

    if (~isempty(currentIndices))
        endTimePoints(i) = max(currentIndices);
    end

    if (trackLengths(i) > trackLengthThreshold || endTimePoints(i) == numFrames)
        trackLabels(i) = currentTrackLabel;
        currentTrackLabel = currentTrackLabel + 1;
    end
end

%% save result images
parfor i=1:numFrames

    %% initialize new result image
    resultImageCleaned = zeros(size(resultImages{i}));

    %% only add objects that were tracked sufficiently long
    validIndices = find(d_orgs(:,i,1) > 0 & (trackLengths > trackLengthThreshold | endTimePoints == numFrames));
    for j=validIndices'
        resultImageCleaned(d_orgs(j,i,3), d_orgs(j,i,4), d_orgs(j,i,5)) = trackLabels(j);
    end

    %% dilate images again for better visibility
    resultImageCleaned = imdilate(resultImageCleaned, strel('sphere', 3));

    %% save the final result images
    saveastiff(uint16(resultImageCleaned), [outputFolderTracked strrep(inputFilesNucleiBasins(i).name, '.tif', '_Tracked.tif')], options);
    saveastiff(uint16(max(resultImageCleaned, [], 3)), [outputFolderTrackedMaxProj strrep(inputFilesNucleiBasins(i).name, '.tif', '_TrackedMaxProj.tif')], options);
end