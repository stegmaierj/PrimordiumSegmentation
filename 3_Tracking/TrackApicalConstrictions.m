%%
% Primordium Segmentation & Tracking.
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

%% specify the parameters
prompt = {'Template radius (px):','Forward Mode (0/1):'};
dlgtitle = 'Please specify the input parameters.';
dims = [1 35];
definput = {'20','0'};
answer = inputdlg(prompt,dlgtitle,dims,definput)

templateRadius = str2double(answer{1});
forwardMode = str2double(answer{2}) > 0;

%% load the input maximum projections (TODO: replace with file open dialog!)
addpath('../ThirdParty/saveastiff_4.3/')
addpath('../ThirdParty/ginputc/')
inputImageXY = loadtiff('/Users/jstegmaier/ScieboDrive/Projects/2022/KnautNYU_PrimordiumCellSegmentation/20211017_Stiching_EGFP Caax-H2a-mCherry_crop_Membranes_MaxProjXY.tif');
inputImageXZ = loadtiff('/Users/jstegmaier/ScieboDrive/Projects/2022/KnautNYU_PrimordiumCellSegmentation/20211017_Stiching_EGFP Caax-H2a-mCherry_crop_Membranes_MaxProjXZ.tif');

%% extract image information
imageSize = size(inputImageXY);
imageSizeXZ = size(inputImageXZ);
numFrames = imageSize(3);

%% specify selection frame
if (forwardMode == true)
    selectionFrame = 1;
else
    selectionFrame = imageSize(3);
end

%% create meshgrid for easier coordinate access
[xpos, ypos, zpos] = meshgrid(1:imageSize(2), 1:imageSize(1), 1:imageSize(3));

%% open annotation figure for landmark selection
figure(1); clf;
colormap gray;
%%set(gcf, 'units','normalized','outerposition',[0 0 1 1], 'Color', [0,0,0]);
colordef white;

subplot(2,1,1); hold on;
imagesc(inputImageXY(:,:,selectionFrame));
axis equal;
axis tight;

subplot(2,1,2); hold on;
imagesc(inputImageXZ(:,:,selectionFrame));
axis equal;
axis tight;

%% initialize the landmark array
apicalConstrictionPoints = [];
answer = questdlg('Add another point?');

%% while loop for adding points to the tracking
while strcmp(answer, 'Yes')

    %% initialize the next point
    newPoint = zeros(1,3);

    %% get point in the XY projection
    subplot(2,1,1);
    %[newPoint(1), newPoint(2)] = ginputc(1, 'Color', 'r', 'LineWidth', 2);
    newPoint(1:2) = ginput(1);

    %[xi, yi] = ginputc(1, 'Color', 'r', 'LineWidth', 2)

    subplot(2,1,1);
    plot(newPoint(1), newPoint(2), '*r');

    %% automatically find the point in the XZ projection
    subplot(2,1,2);
    %xzPoint = ginput(1);
    %newPoint(3) = round(sum((1:size(inputImageXZ, 1)) .* double(inputImageXZ(:, round(newPoint(1)), selectionFrame))') ./ sum(double(inputImageXZ(:, round(newPoint(1)), selectionFrame))));
    [maxValue, maxIndex] = max(imgaussfilt(double(inputImageXZ(:, round(newPoint(1)), selectionFrame)), 3));
    newPoint(3) = maxIndex;
    %newPoint(3) = xzPoint(2);

    subplot(2,1,2);
    plot(newPoint(1), newPoint(3), '*r');

    %% add selected point to the queue
    apicalConstrictionPoints = [apicalConstrictionPoints; round(newPoint)];

    %% ask if another point should be added
    answer = questdlg('Add another point?');
end

%% get the number of points
numInputPoints = size(apicalConstrictionPoints, 1);

%% initialize color map
colors = lines(numInputPoints);

%% initialize the tracking array
acPositions = zeros(numInputPoints, numFrames, 3);
acPositions(:, selectionFrame, :) = apicalConstrictionPoints;

%% specify the frame range and offset depending on the tracking mode
if (forwardMode == true)
    frameRange = 1:(numFrames-1);
    offset = 1;
else
    frameRange = numFrames:(-1):2;
    offset = -1;
end

%% process all frames
for i=frameRange

    %% plot the current figure
    figure(1); clf;
    subplot(2,1,1); 
    imagesc(inputImageXY(:,:,i+offset)); hold on;

    subplot(2,1,2);
    imagesc(inputImageXZ(:,:,i+offset)); hold on;

    %% track the individual points in each frame
    for j=1:numInputPoints

        %% get the previous position
        previousPosition = squeeze(round(acPositions(j, i, :)));

        %% identify the template indices (TODO: boundary check!!!)
        rangeX = (previousPosition(1)-templateRadius):(previousPosition(1)+templateRadius);
        rangeY = (previousPosition(2)-templateRadius):(previousPosition(2)+templateRadius);
        %rangeZ = (previousPosition(3)-templateRadius):(previousPosition(3)+templateRadius);
    
        %% refine the identified position by replacing it with the intensity-weighted centroid
        templatePositionsX = xpos(rangeY, rangeX, 1);
        templatePositionsY = ypos(rangeY, rangeX, 1);
        %templatePositionsZ = zpos(1, rangeX, rangeZ);
        templatePatchXY = double(inputImageXY(rangeY, rangeX, i));
        %templatePatchXZ = double(inputImageXZ(rangeZ, rangeX, i));
        refinedPositionX = sum(templatePositionsX(:) .* templatePatchXY(:)) ./ sum(templatePatchXY(:));
        refinedPositionY = sum(templatePositionsY(:) .* templatePatchXY(:)) ./ sum(templatePatchXY(:));
        %refinedPositionZ = sum(templatePositionsZ(:) .* templatePatchXZ(:)) ./ sum(templatePatchXZ(:));

        %% set the updated position
        previousPosition = round([refinedPositionX, refinedPositionY]);
        acPositions(j, i, 1:2) = previousPosition;

        %% extract the refined region    
        rangeX = (previousPosition(1)-templateRadius):(previousPosition(1)+templateRadius);
        rangeY = (previousPosition(2)-templateRadius):(previousPosition(2)+templateRadius);
        %rangeZ = (previousPosition(3)-templateRadius):(previousPosition(3)+templateRadius);
        templatePatchXY = double(inputImageXY(rangeY, rangeX, i));
        %templatePatchXZ = double(inputImageXZ(rangeZ, rangeX, i));
    
        %% find template in the previous/next frame        
        correlationXY = normxcorr2(templatePatchXY, inputImageXY(:,:,i+offset));
        %correlationXZ = normxcorr2(templatePatchXZ, inputImageXZ(:,:,i+offset));
    
        %% only look for peaks in a small neighborhood
        correlationImageXY = zeros(imageSize(1), imageSize(2));
        correlationImageXY(rangeY, rangeX) = correlationXY(rangeY+templateRadius, rangeX + templateRadius);
        [ypeakXY,xpeakXY] = find(correlationImageXY==max(correlationImageXY(:)));

        %correlationImageXZ = zeros(imageSizeXZ(1), imageSize(2));
        %correlationImageXZ(rangeZ, rangeX) = correlationXZ(rangeZ+templateRadius, rangeX + templateRadius);
        %[ypeakXZ,xpeakXZ] = find(correlationImageXZ==max(correlationImageXZ(:)));
    
        %% save the identified position
        acPositions(j, i+offset, 1:2) = [xpeakXY, ypeakXY];
        %acPositions(j, i+offset, 3) = ypeakXZ;

        %% find the z-position as the maximum intensity spot (TODO: potentially replace with another template matching iteration on the XZ plane)
        [maxValue, maxIndex] = max(imgaussfilt(double(inputImageXZ(:, round(acPositions(j, i+offset, 1)), i+offset)), 3));
        acPositions(j, i+offset, 3) = maxIndex; %round(sum((1:size(inputImageXZ, 1)) .* double(inputImageXZ(:, round(acPositions(j, i-1, 1)), i-1))') ./ sum(double(inputImageXZ(:, round(acPositions(j, i-1, 1)), i-1))));

        %% draw circle at the identified location
        subplot(2,1,1);
        drawcircle(gca, 'Center', squeeze(acPositions(j, i+offset, 1:2))', 'Radius', templateRadius, 'FaceAlpha',0, 'Color', colors(j,:));

        subplot(2,1,2);
        plot(xpeakXY, acPositions(j, i+offset, 3), '*r', 'Color', colors(j,:));
        drawcircle(gca, 'Center', squeeze(acPositions(j, i+offset, [1,3]))', 'Radius', templateRadius, 'FaceAlpha',0, 'Color', colors(j,:));
    end

    drawnow;
end

%% plot the identified trajectories
figure(3); clf; hold on;
for i=1:numPositions
    plot3(1:numFrames, acPositions(i,:,1), acPositions(i,:,2), '-k', 'Color', colors(i,:));
end

%% TODO: export tracks in a proper format to merge/compare them with the cell tracks