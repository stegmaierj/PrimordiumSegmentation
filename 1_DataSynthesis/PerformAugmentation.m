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

function [maskImage, meanImage, stdImage] = PerformAugmentation(maskImage, meanImage, stdImage, settings)
    
    %% perform rotation
    angleX = rand * 360;
    angleY = rand * 360;
    angleZ = rand * 360;
    
    if (settings.rotateX == true)
        maskImage = imrotate3(maskImage, angleX, [1,0,0], 'nearest', 'loose', 'FillValues', 0);
        meanImage = imrotate3(meanImage, angleX, [1,0,0], 'linear', 'loose', 'FillValues', 0);
        stdImage = imrotate3(stdImage, angleX, [1,0,0], 'linear', 'loose', 'FillValues', 0);
    end
    
    if (settings.rotateY == true)
        maskImage = imrotate3(maskImage, angleY, [0,1,0], 'nearest', 'loose', 'FillValues', 0);
        meanImage = imrotate3(meanImage, angleY, [0,1,0], 'linear', 'loose', 'FillValues', 0);
        stdImage = imrotate3(stdImage, angleY, [0,1,0], 'linear', 'loose', 'FillValues', 0);
    end
    
    if (settings.rotateZ == true)
        maskImage = imrotate3(maskImage, angleZ, [0,0,1], 'nearest', 'loose', 'FillValues', 0);
        meanImage = imrotate3(meanImage, angleZ, [0,0,1], 'linear', 'loose', 'FillValues', 0);
        stdImage = imrotate3(stdImage, angleZ, [0,0,1], 'linear', 'loose', 'FillValues', 0);
    end
    
    %% perform scaling
    scaleX = 1 + settings.scaleFactor * (rand - 0.5);
    scaleY = 1 + settings.scaleFactor * (rand - 0.5);
    scaleZ = 1 + settings.scaleFactor * (rand - 0.5);
    
    maskImage = imresize3(maskImage, round(size(maskImage) .* [scaleX, scaleY, scaleZ]), 'nearest');
    meanImage = imresize3(meanImage, round(size(meanImage) .* [scaleX, scaleY, scaleZ]), 'linear');
    stdImage = imresize3(stdImage, round(size(stdImage) .* [scaleX, scaleY, scaleZ]), 'linear');
    
    imageSize = size(maskImage);
    regionProps = regionprops((meanImage+stdImage)>0, 'BoundingBox');
    aabb = regionProps(1).BoundingBox;
    regionPadding = 0;
    rangeX = round(max(1, aabb(1)-regionPadding)):round(min(imageSize(2), aabb(1)+aabb(4)+regionPadding));
    rangeY = round(max(1, aabb(2)-regionPadding)):round(min(imageSize(1), aabb(2)+aabb(5)+regionPadding));
    rangeZ = round(max(1, aabb(3)-regionPadding)):round(min(imageSize(3), aabb(3)+aabb(6)+regionPadding));
    
    maskImage = maskImage(rangeY, rangeX, rangeZ);
    meanImage = meanImage(rangeY, rangeX, rangeZ);
    stdImage = stdImage(rangeY, rangeX, rangeZ);
end