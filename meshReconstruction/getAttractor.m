function attractor = getAttractor(SkipSlice,PointCloudAll,SliceNumAll,attractorPerAxis)
AxisNum = length(SliceNumAll);

for i = 1:AxisNum
    if(~(any(SkipSlice == i)))
    
        if(i >= 8)
            attractorNum(i) = attractorPerAxis*2;
        else
            attractorNum(i) = attractorPerAxis;
        end
    else
        attractorNum(i) = 0;

    end
    
    if(i == 1)
        attractorAll(i) = attractorNum(i);
    else
        attractorAll(i) = attractorAll(i-1) + attractorNum(i);
    end
end

attractor = zeros(attractorAll(length(attractorAll)),3);

 % a total random way to pick up points
 for i = 1:AxisNum
     if(~(any(SkipSlice == i)))
        if(i == 1)
            Point = PointCloudAll(1:SliceNumAll(i),:);
            res = randperm(SliceNumAll(i));
            res = res(1:attractorNum(i));
            attractor(1 : attractorAll(i),: ) = Point(res,:);
        else
            Point = PointCloudAll(SliceNumAll(i-1)+1:SliceNumAll(i),:);
            res = randperm(SliceNumAll(i) - SliceNumAll(i-1));
            res = res(1:attractorNum(i));
            attractor( attractorAll(i-1)+1 : attractorAll(i),: ) = Point(res,:);
        end
     end
end

% attracotr points have a fixed distance from each other
%for i = 1:AxisNum
%    if(i == 1)
%        distance = ceil(SliceNumAll(i) / attractorNum(i));
%        Point = PointCloudAll(1:SliceNumAll(i),:);
%        res = ceil(linspace(1,SliceNumAll(i)-distance,attractorNum(i)));
%        attractor( 1 : attractorAll(i),: ) = Point(res,:);
%    else
%        distance = ceil((SliceNumAll(i) - SliceNumAll(i-1)) / (attractorNum(i)));
%        Point = PointCloudAll(SliceNumAll(i-1)+1:SliceNumAll(i),:);
%        res = ceil(linspace(1,SliceNumAll(i)-SliceNumAll(i-1)-distance,attractorNum(i)));
        
%        attractor( attractorAll(i-1)+1 : attractorAll(i),: ) = Point(res,:);
%    end       
%end