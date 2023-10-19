
areas = [];
for i=1:length(regionPropsTransformed)

    currentAreas = [];
    for j=1:length(regionPropsTransformed{i})   
        currentAreas = [currentAreas; regionPropsTransformed{i}(j).Area];
    end

    areas = [areas; currentAreas];
end

figure;
histogram(areas);

mean(areas)
std(areas)