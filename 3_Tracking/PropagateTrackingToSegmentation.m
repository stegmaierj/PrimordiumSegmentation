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

inputFolderSeg = 'X:\Projects\KnautNYU_PrimordiumCellSegmentation\2022_08_12_NucleiResults_Cellpose+Watershed\20211017_Stiching_EGFP_Caax-H2a-mCherry_crop\';
inputFolderTra = 'C:\Users\stegmaier\Downloads\ManualInitTest\Nuclei_Tracked\';

outputFolder = 'C:\Users\stegmaier\Downloads\ManualInitTest\Nuclei_Tracked_Final\';

inputFilesSeg = dir([inputFolderSeg '*.tif']);
inputFilesTra = dir([inputFolderTra '*.tif']);

clear options;
options.overwrite = true;
options.compress = 'lzw';

parfor i=1:length(inputFilesSeg)

    currentSegImage = loadtiff([inputFolderSeg inputFilesSeg(i).name]);
    currentTraImage = loadtiff([inputFolderTra inputFilesTra(i).name]);

    resultImage = zeros(size(currentTraImage));

    regionPropsSeg = regionprops(currentSegImage, 'Area', 'PixelIdxList');
    regionPropsTra = regionprops(currentTraImage, 'Area', 'PixelIdxList');

    for j=1:length(regionPropsTra)

        if (regionPropsTra(j).Area <= 0)
            continue;
        end
        
        segLabels = unique(currentSegImage(regionPropsTra(j).PixelIdxList));
        segLabels(segLabels == 0) = [];

        if (isempty(segLabels))
            continue;
        end

        for k=segLabels'
            resultImage(regionPropsSeg(k).PixelIdxList) = j;
        end
    end

    saveastiff(uint16(resultImage), [outputFolder inputFilesTra(i).name], options);
end