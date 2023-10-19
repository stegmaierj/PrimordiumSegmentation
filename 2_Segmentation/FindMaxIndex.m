function [matchIndex, matchOverlap] = FindMaxIndex(potentialMatches)

    potentialLabels = unique(potentialMatches);

    potentialLabels(potentialLabels == 0) = [];
    potentialCounts = zeros(size(potentialLabels));

    for k=1:length(potentialLabels)
        potentialCounts(k) = sum(potentialMatches == potentialLabels(k));
    end

    [~, maxIndex] = max(potentialCounts);
    matchIndex = potentialLabels(maxIndex);
    matchOverlap = potentialCounts(maxIndex);
end