function distance = DistanceBetweenPoints(p1,p2)

    if(size(p2,1))==1
        distance = sqrt(sum((p1 - repmat(p2,size(p1,1),1)) .^2,2));
    else
        distance = sqrt(sum((p2 - repmat(p1,size(p2,1),1)) .^2 ,2));
        %   distance = sqrt((p1(1) - rep2(1))^2 + (p1(2)-p2(2))^2 + (p1(3)-p2(3))^2);
    end
end