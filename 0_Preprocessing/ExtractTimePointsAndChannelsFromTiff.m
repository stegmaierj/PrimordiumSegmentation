
inputFileName = '/Users/jstegmaier/Downloads/20231026_prim-mCherry_H2A-GFP_tracking_0.5 hour.tif';
%inputFileName = '/Users/jstegmaier/Downloads/20231025_EGFP-CaaX_H2A-mCherry_tracking_crop_0.5hour.tif';

inputImage = loadtiff(inputFileName);

[folder, file, ext] = fileparts(inputFileName);

outputFolderNuclei = ['/Users/jstegmaier/Downloads/' file '_Nuclei/'];
outputFolderMembranes = ['/Users/jstegmaier/Downloads/' file '_Membranes/'];
if (~exist(outputFolderNuclei, 'dir')); mkdir(outputFolderNuclei); end
if (~exist(outputFolderMembranes, 'dir')); mkdir(outputFolderMembranes); end

numChannels = 2;
numSlices = 70;
numFrames = 20;
imageHeight = 173;
imageWidth = 1262;

currentImageCh1 = zeros(imageHeight, imageWidth, numSlices);
currentImageCh2 = zeros(imageHeight, imageWidth, numSlices);

clear options;
options.overwrite = true;
options.compress = 'lzw';

for i=1:numFrames
    currentImageCh1 = zeros(imageHeight, imageWidth, numSlices);
    currentImageCh2 = zeros(imageHeight, imageWidth, numSlices);

    for j=1:2:(numSlices*numChannels)
        imageIndexCh1 = (i-1)*numSlices*numChannels + j;
        imageIndexCh2 = (i-1)*numSlices*numChannels + j + 1;

        currentImageCh1(:,:,(j+1) / 2) = inputImage(:,:,imageIndexCh1);
        currentImageCh2(:,:, (j+1) / 2) = inputImage(:,:,imageIndexCh2);

    end

    saveastiff(uint16(currentImageCh2), sprintf('%s/%s_c=%02d_t=%04d.tif', outputFolderNuclei, file, 0, i), options);
    saveastiff(uint16(currentImageCh1), sprintf('%s/%s_c=%02d_t=%04d.tif', outputFolderMembranes, file, 1, i), options);
    
    test = 1;
end


figure(1);
subplot(1,2,1)
imagesc(inputImage(:,:,31));

subplot(1,2,2)
imagesc(inputImage(:,:,32));