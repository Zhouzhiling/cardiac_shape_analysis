function [df, diffVector] = calculateDeformFactor(mesh,attractor)
    pointNum = length(mesh.vertices(:,1));
    diffVector = zeros(pointNum,3);
    df = zeros(1,pointNum);
    attractorNum = length(attractor(:,1));
    for i = 1:pointNum
        curPoint = mesh.vertices(i,:);
        [distance, closePoint] = findClosestPoint(curPoint,attractor,attractorNum);
        diffVector(i,:) = closePoint - curPoint;
        %dir = diffVector * mesh.normal(i,:)';
        
        %if(DistanceBetweenPoints(centerPoint,curPoint) > DistanceBetweenPoints(centerPoint,closePoint))
        %if(dir < 0)
        %    df(i) = -distance;
        %else
        %    df(i) = distance;
        %end
    end
    %diffVector = rescale(diffVector)
    %width = length(df(1,:));
    %height = length(df(:,1));
    %df = reshape(df,1,width*height);
    %df = mapminmax(df,0,1);
    %df = reshape(df,height,width);
   
end

function res = rescale(input)
    res = zeros(length(input(:,1)),length(input(1,:)));
    num = length(input(:,1));
    for i = 1:num
        scale = sqrt(input(i,1).^2 + input(i,2).^2 + input(i,3).^2 );
        res(i,:) = input(i,:) ./ scale;
    end
end

function [minD, cp] = findClosestPoint(cup,attractor,attractorNum)
    minD = 5000;
    for j = 1:attractorNum
        dist = DistanceBetweenPoints(attractor(j,:),cup);
        if dist < minD
            minD = dist;
            cp = attractor(j,:);
        end 
    end    
end

function distance = DistanceBetweenPoints(p1,p2)
    distance = sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (p1(3)-p2(3))^2);
end