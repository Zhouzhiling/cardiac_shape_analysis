function [pos] = register(rl,ap,fh,ImageOrientation,ImagePosition,Resolution)
%---------------------------------------
%Convert from RL,AP,FH coordinates to Segment internal coordinate system.

  Resolution = Resolution(1);
  xdir = ImageOrientation(4:6)';
  ydir = ImageOrientation(1:3)';  
  zdir = cross(...
    ImageOrientation(1:3),...
    ImageOrientation(4:6))';    
  
  rl = rl(:)'-ImagePosition(1); %Translate corner of box to 0,0,0
  ap = ap(:)'-ImagePosition(2);
  fh = fh(:)'-ImagePosition(3);
  
  pos = [xdir ydir -zdir]\[rl;ap;fh];
%  pos(1,:) = pos(1,:)+1;
%  pos(2,:) = pos(2,:)+1;
%  pos(3,:) = pos(3,:)+1;  