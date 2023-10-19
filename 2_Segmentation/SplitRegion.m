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

function splitRegion = SplitRegion(imageRegion, currAABB)



    seedCandidates = nextLabels(validMatchesNext2Curr);

    labelImage = zeros(size(imageRegion));
    for k=seedCandidates'

        seedImage = (double(nextSegmentation(rangeX, rangeY, rangeZ)) .* (double(currSegmentation(rangeX, rangeY, rangeZ)) == currentId)) == k;
        [xpos, ypos, zpos] = ind2sub(size(seedImage), find(seedImage(:) > 0));

        labelImage(round(mean(xpos)), round(mean(ypos)), round(mean(zpos))) = 1;
    end
    imageRegion = imimposemin(imageRegion, labelImage);

    currenRegionWs = watershed(imageRegion);
    splitImage = double(currenRegionWs > 0) .* double(currSegmentation(rangeX, rangeY, rangeZ) == currentId);
    imagesc(splitImage(:,:,round(size(splitImage,3)/2)));

    %% replace the original segment with the split one
    %% TODO
    oldRegionContent = currSegmentation(rangeX, rangeY, rangeZ);
    oldRegionContent(oldRegionContent == currentId) = 0;

    localRegionProps = regionprops(splitImage, 'PixelIdxList');
    for k=1:length(localRegionProps)
        if (k==1)
            oldRegionContent(localRegionProps(k).PixelIdxList) = currentId;
        else
            oldRegionContent(localRegionProps(k).PixelIdxList) = nextNewLabel;
            nextNewLabel = nextNewLabel+1;
        end
    end

    currSegmentationCorrected(rangeX, rangeY, rangeZ) = oldRegionContent;

end