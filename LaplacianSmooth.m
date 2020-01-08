function dataAfterSmooth = LaplacianSmooth(data, LaplacianSmoothTimes)

% if nargin < 2 
%    LaplacianSmoothTimes = 1; 
% end
% my own function to smooth the point cloud

% this assignment to data is for test only 
% data = [xi', yi', zi'./10];
dataAfterSmooth = data;
pointNumber = length(data(:,1));
sliceNumber = length(unique(data(:,3)));

% 求每一层的元素个数
pointNumberPerSlice = zeros(1,sliceNumber);
pointNumberPerSingleSlice = zeros(1,sliceNumber);
% pointNumberPerSlice = [];
% pointNumberPerSingleSlice = [];

curData = [];
preData = [];
prepreData = [];

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

% 遍历每一层（每一个z值
% for h = 1:sliceNumber

for k = 1:LaplacianSmoothTimes
    for h = 1:sliceNumber
        if(h == 1)
            % 1
            % pointNumberPerSlice(h)
            prepreData = data(1:pointNumberPerSlice(h),:);
            preData = data(1:pointNumberPerSlice(h),:);
            curData = data(1:pointNumberPerSlice(h),:);
            % avoid the first layer from laplacian
        else
            prepreData = preData;
            preData = curData;
            curData = data(pointNumberPerSlice(h-1)+1:pointNumberPerSlice(h),:);
        end

        if(h > 2)
            % start smoothing the preData slice
            toUpdatedSlice = preData;
            % traverse all points in this slice get the new numbers
            for i = 1:pointNumberPerSingleSlice(h-1)
                preData(i,:) = findClosetPoint(preData(i,:),preData, prepreData,curData);
                % preData(i,:) = findClosetPoint(preData(i,:),preData, prepreData,curData);
            end

            % write data back
            dataAfterSmooth(pointNumberPerSlice(h-2)+1:pointNumberPerSlice(h-1),:) = preData;
            data(pointNumberPerSlice(h-2)+1:pointNumberPerSlice(h-1),:) = preData;
        end
    end
end
%dataAfterSmooth(:,3) = dataAfterSmooth(:,3).*10;

%attractorPoints = randomlyPickPoints(data, pointNumberPerSingleSlice, pointNumberPerSlice, sliceNumber, 4);


    