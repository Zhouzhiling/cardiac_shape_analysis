function [pos] = xyz2rlapfh(x,y,z,ImageOrientation,ImagePosition,Resolution,SliceThickness)
%------------------------------------
%Converts from segment coordinate system to RL,AP,FH coordinate system.
% for test
    % x = xi(i);
    % y = yi(i)
    % z = zi(i);
 
  Resolution = Resolution(1);

  zdir = cross(...
    ImageOrientation(1:3),...
    ImageOrientation(4:6)); 

  x = (x(:)-1)*Resolution;
  y = (y(:)-1)*Resolution;
  z = (z(:)-1)*SliceThickness;
  
  pos = repmat(ImagePosition,length(x),1)+...
    repmat(ImageOrientation(4:6),length(x),1).*repmat(x,1,3)+...
    repmat(ImageOrientation(1:3),length(x),1).*repmat(y,1,3)-... %was minus
    repmat(zdir,length(x),1).*repmat(z,1,3);    
