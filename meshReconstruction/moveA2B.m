function result = moveA2B(A,B)


%A = LAX2;
%B = PointCloud_o;
%A = vertices';
%B = EpicardiumBeforeSmooth;
% result = [];

meanAX = mean(A(:,1));
meanAY = mean(A(:,2));
meanAZ = mean(A(:,3));

meanBX = mean(B(:,1));
meanBY = mean(B(:,2));
meanBZ = mean(B(:,3));

moveX = mean(B(:,1)) - mean(A(:,1));
moveY = mean(B(:,2)) - mean(A(:,2));
moveZ = mean(B(:,3)) - mean(A(:,3));

result = [A(:,1) + moveX, A(:,2) + moveY, A(:,3) + moveZ];


