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

function [position] = FindValidSamplingLocation(maskImage, minDistance)

    %% create distance map of current mask
    distanceMapImage = bwdist(maskImage > 0);
    
    %% idenfify the valid positions to sample from
    validLocations = find(distanceMapImage(:) > minDistance);

    %% draw one random position for placing the new object at
    position = zeros(3,1);
    if (~isempty(validLocations))
        [position(1), position(2), position(3)] = ind2sub(size(maskImage), validLocations(randperm(length(validLocations),1)));
    else
        position = [];
        disp('No remaining position available for current min distance!');
    end
end