function dataAfter = addUpandDownLimit(data)

sliceHeight = 8;
% data = EndocardiumAfterSmooth;
sliceNumber = length(unique(data(:,3)));

% 求每一层的元素个数
pointNumberPerSlice = zeros(1,sliceNumber);
pointNumberPerSingleSlice = zeros(1,sliceNumber);

curHeight = data(1,3);
curNumber = length(find(data(:,3)==curHeight));
pointNumberPerSlice(1) = curNumber;
pointNumberPerSingleSlice(1) = curNumber;

for h = 2:sliceNumber
    firstSlice = pointNumberPerSlice(h-1) + 1;
    curHeight = data(firstSlice,3);
    curNumber = length(find(data(:,3)==curHeight));
    addedNumber = curNumber + pointNumberPerSlice(h-1);
    pointNumberPerSlice(h) = addedNumber;
    pointNumberPerSingleSlice(h) = curNumber;
end

totalX = 0;
totalY = 0;

% calculate the lower bound
for i = 1:pointNumberPerSingleSlice(1)
    totalX = totalX + data(i, 1);
    totalY = totalY + data(i, 2);
end

averageX = totalX / pointNumberPerSingleSlice(1);
averageY = totalY / pointNumberPerSingleSlice(1);
averageZ = data(1,3);
lowerBound = [averageX, averageY, averageZ - sliceHeight];

% calculate the upper bound
totalX = 0;
totalY = 0;

for i = pointNumberPerSlice(sliceNumber-1)+1:pointNumberPerSlice(sliceNumber)
    totalX = totalX + data(i, 1);
    totalY = totalY + data(i, 2);
end
averageX = totalX / pointNumberPerSingleSlice(sliceNumber);
averageY = totalY / pointNumberPerSingleSlice(sliceNumber);
averageZ = 2.*data(i,3) - data(pointNumberPerSlice(sliceNumber-1),3);
upperBound = [averageX, averageY, averageZ];
dataAfter = [lowerBound;data;upperBound];






