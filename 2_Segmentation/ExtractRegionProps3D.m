
addpath('../ThirdParty')
addpath('../ThirdParty/saveastiff_4.3');

outputFolder = 'X:\Projects\KnautNYU_PrimordiumCellSegmentation\2021_10_08\20210506_Stiching_EGFP_Caax-H2a-mCherry\Membrane_RegionProps\';
if (~isfolder(outputFolder)); mkdir(outputFolder); end

inputFolder = 'X:\Projects\KnautNYU_PrimordiumCellSegmentation\2021_10_08\20210506_Stiching_EGFP_Caax-H2a-mCherry\Membranes_Tracked\';
inputFiles = dir([inputFolder '*.tif']);

zspacing = 2.4615; %% 0.325 x 0.325 x 0.8 µm 

specifiers = 'id;volume;xpos;ypos;zpos;equivDiameter;princialAxisLength1;principalAxisLength2;principalAxisLength3;orientation1;orientation2;orientation3;convexVolume;solidity;surfaceArea';

parfor i=1:length(inputFiles)
    
    currentImage = loadtiff([inputFolder inputFiles(i).name]);
    
    currentImage = imresize3(currentImage, [size(currentImage,1), size(currentImage,2), round(zspacing * size(currentImage,3))], 'nearest');
    
    currentRegionProps = regionprops3(currentImage, 'Centroid', 'PrincipalAxisLength', 'Volume', 'SurfaceArea', 'ConvexVolume', 'EquivDiameter', 'Solidity', 'Solidity', 'Orientation');
    
    resultTable = table2array(currentRegionProps);
    
    volumeIndex = 1;
    centroidIndex = 2:4;
    equivDiameterIndex = 5;
    principalAxisLength = 6:8;
    orientationIndex = 9:11;
    convexVolumeIndex = 12;
    solidityIndex = 13;
    surfaceAreaIndex = 14;
    
    ids = 1:size(resultTable,1);
    
    resultTable = [ids', resultTable];
    resultTable(resultTable(:,2) == 0, :) = [];
    
    outputFileName = [outputFolder strrep(inputFiles(i).name, '.tif', '_RegionProps.csv')];
    dlmwrite(outputFileName, resultTable, ';');
    prepend2file(specifiers, outputFileName, true);
    test = 1;
end