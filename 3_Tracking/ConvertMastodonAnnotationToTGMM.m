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

%% NOTE: This is just a correction factor here, as the initial image in Fiji had this weird spacing contained.
%% So this is not the real physical spacing but is used to convert the coordinates to image coordinates
pixelSpacing = [0.039370079301150995 0.039370079301150995 0.039370079301150995];

%% setup input paths
addpath('../ThirdParty')
addpath('../ThirdParty/saveastiff_4.3');
addpath('../ThirdParty/matlab-mastodon-importer/src');

%% specify the output directory
outputRoot = uigetdir('C:\Users\stegmaier\Downloads\GlobalOutputTest\', 'Please select the output root directory ...');
outputRoot = [outputRoot filesep];
if (~isfolder(outputRoot)); mkdir(outputRoot); end

%% specify the input directory
inputFolderTrackedNuclei = [outputRoot 'Nuclei_Tracked' filesep];
if (~isfolder(inputFolderTrackedNuclei))
    inputFolderTrackedNuclei = uigetdir('X:\Projects\KnautNYU_PrimordiumCellSegmentation\2021_10_08\20210506_Stiching_EGFP_Caax-H2a-mCherry\Nuclei_Tracked\', 'Please select the input directory (should contain tracked 3D tiffs for each frame)');
    inputFolderTrackedNuclei = [inputFolderTrackedNuclei filesep];
else
   disp(['Using tracked nuclei located in ' inputFolderTrackedNuclei]); 
end
inputFiles = dir([inputFolderTrackedNuclei '*.tif']);

%% specify the output directory
outputFolderTracked = [outputRoot 'Nuclei_Tracked_Corrected' filesep];
if (~isfolder(outputFolderTracked)); mkdir(outputFolderTracked); end

%% specify the output directory
outputFolderTrackedMaxProj = [outputRoot 'Nuclei_Tracked_Corrected_MaxProj' filesep];
if (~isfolder(outputFolderTrackedMaxProj)); mkdir(outputFolderTrackedMaxProj); end

%% specify the output directory
outputFolderTGMM = [outputRoot 'TGMM_Corrected' filesep];
if (~isfolder(outputFolderTGMM)); mkdir(outputFolderTGMM); end

%% specify the input mastodon file
[inputFileMastodon, inputPathMastodon]  = uigetfile('*.mastodon', 'Select the mastodon file you want to apply to the selected images');

%% bring mastodon file into a d_orgs like format for easier application to the image data
[G, metadata, tss] = import_mastodon( [inputPathMastodon inputFileMastodon] );

nodes = table2array(G.Nodes(:,1:5));
nodes(:,1) = nodes(:,1)+1;
nodes(:,5) = nodes(:,5)+1;

edges = table2array(G.Edges);

%% sort the edges based on their start time point
%% otherwise it can happen that objects are added later and start a new track
edges(:,end+1) = 0;
numDivisionsExpected = 0;
for i=1:size(edges,1)
    
    currentEdge = edges(i,:);
    
    object1 = nodes(currentEdge(1), :);

    %% check if two daughters exist
    numSuccessors = find(edges(:,1) == currentEdge(1));

    if (length(numSuccessors) > 1)
        disp('Division Detected!');
        numDivisionsExpected = numDivisionsExpected + 1;
    end
   
    timePoint1 = object1(5);
    edges(i,end) = timePoint1;
end

[timePoints, edgeOrder] = sort(edges(:,end), 'ascend');
edges = edges(edgeOrder, :);
numTimePoints = max(timePoints) + 1;

visitedNodes = zeros(size(nodes,1),1);
visitedEdges = zeros(size(edges,1),1);

%% fill d_orgs structure with all nodes (will be very sparse)
d_orgs = zeros(size(nodes,1), numTimePoints, 7);
for i=1:size(nodes,1)
    currentTimePoint = nodes(i,5);
    d_orgs(i, currentTimePoint, 1) = nodes(i,1);
    d_orgs(i, currentTimePoint, 2) = 1;
    d_orgs(i, currentTimePoint, 3:5) = round(nodes(i,2:4) ./ pixelSpacing);
end

%% put corresponding cells to the same row and handle cell divisions
numDivisionsHandled = 0;
deletionIndices = zeros(size(d_orgs,1),1);
for i=1:size(edges,1)

    %% get the current edge
    currentEdge = edges(i,:);
    
    %% get the object ids linked by the edge
    object1 = nodes(currentEdge(1), :);
    object2 = nodes(currentEdge(2), :);
    
    %% get the time points of the two objects
    timePoint1 = object1(5);
    timePoint2 = object2(5);
    
    %% get the index of the objects in d_orgs
    object1Index = find(d_orgs(:,timePoint1,1) == currentEdge(1));
    object2Index = find(d_orgs(:,timePoint2,1) == currentEdge(2));

    %% check how many edges are affected by a change of object 2
    affectedEdgesCurrentFrame = find(edges(:,1) == object1Index & ~visitedEdges);
    affectedEdgesNextFrame = find(edges(:,1) == object2Index & ~visitedEdges);

    currentLineageId = d_orgs(object1Index, timePoint1, 7);
    if (currentLineageId == 0)
        currentLineageId = max(d_orgs(:, :, 7), [], 'all') + 1;
        d_orgs(object1Index, timePoint1, 7) = currentLineageId;
    end
    
    %% handle movement and division separately
    if (length(affectedEdgesCurrentFrame) == 1)

        %% copy the successor to the same row in d_orgs and delete the second entry
        d_orgs(object1Index, timePoint2, :) = d_orgs(object2Index, timePoint2, :);
        d_orgs(object2Index, timePoint2, :) = 0;
        edges(affectedEdgesNextFrame, 1) = object1Index;

        %% set the predecessor id
        d_orgs(object1Index, timePoint2, 6) = object1Index;
        d_orgs(object1Index, timePoint2, 7) = currentLineageId;

        deletionIndices(object2Index) = 1;

    %% if there are two edges present in the current frame that link to object 1, a division is present.
    %% in the case of divisions, the two new IDs are used instead of copying it to the row of the predecessor
    else

        %% get the edges that belong to the division
        edge1 = edges(affectedEdgesCurrentFrame(1),:);
        edge2 = edges(affectedEdgesCurrentFrame(2),:);

        %% set the predecessor index
        d_orgs(edge1(2), timePoint2, 6) = object1Index;
        d_orgs(edge2(2), timePoint2, 6) = object1Index;

        %% set the lineage id
        d_orgs(edge1(2), timePoint2, 7) = currentLineageId;
        d_orgs(edge2(2), timePoint2, 7) = currentLineageId;

        %% sanity check variable to see if all divisions were handled
        numDivisionsHandled = numDivisionsHandled + 1;
    end

    visitedEdges(i) = 1;
end

%% remove empty rows
d_orgs(deletionIndices > 0, :, :) = 0;

%% figure for debugging
% figure(1);
% for i=1:numTimePoints
%     validIndices = find(d_orgs(:,i,1) > 0); 
%     scatter3(d_orgs(validIndices, i, 4), d_orgs(validIndices, i, 3), d_orgs(validIndices, i, 5), 10, d_orgs(validIndices, i, 7), 'filled');
%     axis equal;
%     view(90, 90);
%     drawnow;
% end

%% renumber indices that belong to a single track
%% all objects in one row have the same id
for i=1:size(d_orgs,1)
    validIndices = squeeze(d_orgs(i,:,1)) > 0;
    d_orgs(i,validIndices,1) = i;
end

%% get the unique labels that are occupied by any objects
%% this is used as a mapping, to use a small number of effective labels
uniqueLabels = unique(d_orgs(:,:,1));

%% define z-spacing (only needed if isotropic images should be processed).
zspacing = 1; %% 0.325 x 0.325 x 0.8 µm

%% specify indices of the features in the region props table
volumeIndex = 1;
centroidIndex = 2:4;
equivDiameterIndex = 5;
principalAxisLength = 6:8;
orientationIndex = 9:11;
convexVolumeIndex = 12;
solidityIndex = 13;
surfaceAreaIndex = 14;

%% initialize the previous labels array
previousLabels = [];

sampleImage = loadtiff([inputFolderTrackedNuclei inputFiles(1).name]);
imageSize = size(sampleImage);

%% process all image files sequentially
for i=1:size(d_orgs,2)
    
    %% initialize an empty image
    resultImage = zeros(imageSize);

    %% add the detections to the result image
    %% instead of the original labels, the mapping to the row number in uniqueLabels 
    %% is used to obtain a smaller set of labels and avoid unoccupied labels.
    validIndices = find(d_orgs(:,i,1) > 0);
    for j=validIndices'
        currentCentroid = squeeze(d_orgs(j,i,3:5));
        resultImage(currentCentroid(2), currentCentroid(1), currentCentroid(3)) = find(uniqueLabels == j);
    end

    %% define output filenames and write the current result file
    outputFileName = [outputFolderTracked strrep(inputFiles(i).name, '.tif', '_ManuallyTracked.tif')];
    outputFileNameMaxProj = [outputFolderTrackedMaxProj strrep(inputFiles(i).name, '.tif', '_ManuallyTracked.tif')];

    clear options;
    options.compress = 'lzw';
    options.overwrite = true;
    saveastiff(uint16(resultImage), outputFileName, options);
    saveastiff(uint16(max(imdilate(resultImage, strel('sphere', 3)), [], 3)), outputFileNameMaxProj, options);

    %% potentially resize the input image in z if needed
    if (zspacing ~= 1)
        resultImage = imresize3(resultImage, [size(resultImage,1), size(resultImage,2), round(zspacing * size(resultImage,3))], 'nearest');
    end

    %% extract the region properties
    currentRegionProps = regionprops3(resultImage, 'Centroid', 'PrincipalAxisLength', 'Volume', 'SurfaceArea', 'ConvexVolume', 'EquivDiameter', 'Solidity', 'Solidity', 'Orientation', 'VoxelIdxList');

    voxelIdxLists = currentRegionProps.VoxelIdxList;
    currentRegionProps.VoxelIdxList = [];

    %% convert region props table to a regular array
    resultTable = table2array(currentRegionProps);

    %% create ids for the current list of objects
    ids = 1:size(resultTable,1);

    %% assemble the results table
    resultTable = [ids', resultTable(:,volumeIndex), round(resultTable(:,centroidIndex)), (i-1)*ones(size(ids')), 100*ones(size(ids')), ids', ids'];

    invalidIndices = find(resultTable(:,2) == 0);
    resultTable(invalidIndices, :) = [];
    voxelIdxLists(invalidIndices) = [];

    %% create TGMM style output xml files for each of the frames
    docNode = com.mathworks.xml.XMLUtils.createDocument('document'); %#ok<JAPIMATHWORKS>
    toc = docNode.getDocumentElement;

    %% add each detection as a fake Gaussian mixture model
    %% only the centroid and IDs are reasonable values, 
    %% the rest are just dummy values from one of the TGMM examples.

    validIndices = find(d_orgs(:,i, 1) > 0);
    for j=validIndices'
        currentGM = docNode.createElement('GaussianMixtureModel');
        currentGM.setAttribute('id', num2str( find(uniqueLabels == j) ));
        currentGM.setAttribute('lineage', num2str( find(uniqueLabels == d_orgs(j, i, 7)) ) );

        %% this was used before to estimate the covariance matrix from the image data.
        %% as we have only point locations here, the covariance is just initialized with dummy values
        %[xpos, ypos, zpos] = ind2sub(size(resultImage), voxelIdxLists{j});
        %W = cov([zpos, ypos, xpos]);
        %WString = strrep(strrep(strrep(strrep(strrep(strrep(num2str(W(:)'), '  ', ' '), '  ', ' '), '  ', ' '), '  ', ' '), '  ', ' '), '  ', ' ');

        %% set the parent id if the predecessor is non-zero
        if (d_orgs(j, i, 6) > 0)
            currentGM.setAttribute('parent', num2str( find(uniqueLabels == d_orgs(j, i, 6)) ));
        else
            currentGM.setAttribute('parent', '-1');
        end
        currentGM.setAttribute('dims', '3');
        currentGM.setAttribute('splitScore', '3');
        currentGM.setAttribute('scale', '1 1 1');
        currentGM.setAttribute('nu', '126.937');
        currentGM.setAttribute('beta', '126.937');

        %% set the mean vector
        currentGM.setAttribute('m', num2str( squeeze(d_orgs(j, i, 3:5))' ));

        %% set the (dummy) covariance matrix
        %if (rank(W) < 3)
            currentGM.setAttribute('W', '1 0 0 0 1 0 0 0 1');
        %else
        %    currentGM.setAttribute('W', WString);
        %end

        %% another bunch of dummy parameters that are from the 
        %% TGMM tracking algorithm but most likely not used by Mastodon
        currentGM.setAttribute('nuPrior', '4');
        currentGM.setAttribute('betaPrior', '0.075017');
        currentGM.setAttribute('alphaPrior', '0');
        currentGM.setAttribute('distMRFPrior', '0');
        currentGM.setAttribute('mPrior', num2str( squeeze(d_orgs(j, i, 3:5))' ));
        currentGM.setAttribute('WPrior', '1 0 0 0 1 0 0 0 1');
        currentGM.setAttribute('svIdx', '0');
        toc.appendChild(currentGM);
    end

    %% write the current result file
    outputFileName = [outputFolderTGMM strrep(inputFiles(i).name, '.tif', '_RegionProps.xml')];
    xmlwrite(outputFileName,docNode);

    %% store the previous labels to see if there's a predecessor available    
    previousLabels = resultTable(:,1);
end
