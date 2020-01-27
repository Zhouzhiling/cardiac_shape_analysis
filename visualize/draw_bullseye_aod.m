% draw bullseye
%-----------------
% Author: zhiling zhou
% Date:   15th Sep, 2018  
%-----------------
% clear

for idx = 1:8
	mean_dif = value(idx,:);
	name = NAME{idx+8};
	rr = drawbullseyeV2(mean_dif,name);
	
	outpath2 = [outroot,name,'_bullseye.png'];
	saveas(gcf,outpath2);
end

%%
close;
figure(1347)
c = createBullseye([0 0.5 1 0; 0.5 1 4 45; 1 1.5 6 0; 1.5 2 6 0]);
set(c,'Color','k','LineWidth',2)

% Example 1 of filling the bullseye, vector by vector
% 16 13 14 15
%fillBullseye([-2],0.5,1,-45,315);
fillBullseye([-1 -4],1,1.5,0,360);
%fillBullseye([2],1.5,2,0,360);
%fillBullseye( ( 1:10:360 )/360,1.5,2,0,360);
uistack(c,'top');