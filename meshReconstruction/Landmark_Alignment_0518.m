%% start from here!!!
clear all;
caseName = 'SC-HF-I-11';
Corresponding = 'IM-0047';

txtPath = strcat('./data/dataset/',caseName,'/contours-manual/IRCCI-expert/');
outputPath = './output/0422_artifical_point_cloud/new524/';
%dcmPath = strcat('./data/alldata/',caseName,'/dcm/');
dcmPath = strcat('F:/Harvard/dcm/',caseName,'/');

% V2.0 more data needed
% SC-HF-I-01    0004     
% SC-HF-I-02    0106
% SC-HF-I-04     116
% SC-HF-I-05       0156
% SC-HF-I-06     180
% SC-HF-I-07     209
% SC-HF-I-08     226
% SC-HF-I-09     241
% SC-HF-I-10     0024
% SC-HF-I-11     0047
% SC-HF-I-12     0062
% SC-HF-I-40     0134

% SC-HF-NI-03    0379
% SC-HF-NI-04    0501
% SC-HF-NI-07    0523
% SC-HF-NI-11   0270
% SC-HF-NI-12   0286
% SC-HF-NI-13   0304
% SC-HF-NI-14   0331
% SC-HF-NI-15   0359
% SC-HF-NI-33   424      
% SC-HF-NI-34   446   
% SC-HF-NI-31   401   
% SC-HF-NI-36   0474


% HYP-01        550    
% HYP-03        0650
% SC-HYP-06     0767
% SC-HYP-07     0007
% SC-HYP-08     0796
% SC-HYP-09     0003
% SC-HYP-10     579
% SC-HYP-11     601
% SC-HYP-12     629
% SC-HYP-37     702
% SC-HYP-38     0734
% SC-HYP-40     755

% SC-N-02     898
% SC-N-03         915
% SC-N-05     963
% SC-N-06     984
% SC-N-07     1009
% SC-N-09     1031
% SC-N-10     851
% SC-N-11     878
% SC-N-40     0944

% read txt file
flagi = 0;
flago = 0;
flagp = 0; 

listing = dir(txtPath);
numofFile = length(listing);
PointCloud = [];
SliceNum = [];
SliceNumTotal = [];
LandMarkE = [];
LandMarkS = [];
% start from the third file, skip the first two files, '.' and '..'
% 1~20      first frame
% 21~40     second frame
% k         ceil(k/20) frame
numO = 0;

for i = 3:numofFile
% for i = 4
    LandExist = false;
    newSliceGet = false;
    fileName = listing(i).name;
    ct = load(strcat(txtPath,fileName));
    IDString = fileName(9:12);
    ID = str2double(fileName(10:12));
    numSlice = ceil(ID/20);
    fileType = fileName(14);
    
    AnchorPath = strcat('./data/SAX_done/',caseName,'/');
    AnchorF = strcat(AnchorPath,'ED-',IDString,'-',caseName,'_anchor.txt');
    
    if(fileType == 'i' && ~(flagi == numSlice))
%        % EndoES
%        numO = numO + 1;
%        SliceNum(numO) = length(ct(:,1));
        type = 1;newSliceGet = true;flagi = numSlice;
%        c = linspace(1,1,length(ct(:,1)));
%        xo = ct(:,1)';
%        yo = ct(:,2)';
%        zo = c';
%        fprintf('Endo - ES - ID is %s\n',IDString);       
        
    elseif(fileType == 'i' && (flagi == numSlice))
        % Endo ED
%        numO = numO + 1;
%        SliceNum(numO) = length(ct(:,1));
        type = 2;newSliceGet = true;
%        c = linspace(1,1,length(ct(:,1)));
%        xo = ct(:,1)';
%        yo = ct(:,2)';
%        zo = c';
%        fprintf('Endo - ED - ID is %s\n',IDString);

    % else if the points represent epicardium  
    elseif(fileType == 'o' && ~(flago == numSlice))
        numO = numO + 1;
        SliceNum(numO) = length(ct(:,1));
        type = 3;newSliceGet = true;flago = numSlice;
        c = linspace(1,1,length(ct(:,1)));
        xo = ct(:,1)';
        yo = ct(:,2)';
        zo = c';
        fprintf('Epi - ED - ID is %s\n',IDString);
        
         if(exist(AnchorF,'file'))
            LandExist = true;
            LM = load(AnchorF);
            LandS = [LM(1,:), zo(1)];
            LandE = [LM(2,:), zo(1)];
            fprintf('Landmark is %s\n',IDString);
        else
            disp('NO');
         end
    else
        newSliceGet = false;
    end
    
  
    % read corresponding dcm
    if(newSliceGet)
        dcmName = strcat(Corresponding,'-',IDString,'.dcm');
        DCM = dicominfo(strcat(dcmPath,dcmName));
        % get information from the dcm file
        ImageOrientation        =    DCM.ImageOrientationPatient';
        PixelSpacing            =    DCM.PixelSpacing;
        SliceThickness          =    DCM.SliceThickness;
        ImagePosition           =    DCM.ImagePositionPatient';
        
        % transform this txt file, keep the position
      %  if(type == 1)
      %      for idx = 1 : length(xi)
      %          newPos = xyz2rlapfh(xi(idx),yi(idx),zi(idx),ImageOrientation,ImagePosition,PixelSpacing,SliceThickness);
      %          PointCloud = [PointCloud; newPos];
      %      end
            
        %elseif(type == 3)
        if(type == 3)
            % epi es
            for idx = 1 : length(xo)
                newPos = xyz2rlapfh(xo(idx),yo(idx),zo(idx),ImageOrientation,ImagePosition,PixelSpacing,SliceThickness);
                PointCloud = [PointCloud; newPos];
            end
            if(LandExist)
                LandSNew = xyz2rlapfh(LandS(1),LandS(2),LandS(3),ImageOrientation,ImagePosition,PixelSpacing,SliceThickness);
                LandENew = xyz2rlapfh(LandE(1),LandE(2),LandE(3),ImageOrientation,ImagePosition,PixelSpacing,SliceThickness);

                LandMarkS = [LandMarkS; LandSNew];
                LandMarkE = [LandMarkE; LandENew];
            end
     %   elseif(type == 2)
     %       % endo ed
     %       for idx = 1 : length(xied)
     %           newPos = xyz2rlapfh(xied(idx),yied(idx),zied(idx),ImageOrientation,ImagePosition,PixelSpacing,SliceThickness);
     % %          PointCloud = [PointCloud; newPos];
     %       end
        end
    end
end


%figure(546)
%scatter3(PointCloud(:,1),PointCloud(:,2),PointCloud(:,3));
% register



for i = 1:length(SliceNum)
    if(i==1)
        SliceNumAll(i) =  SliceNum(i);
    else
        SliceNumAll(i) =  SliceNum(i) + SliceNumAll(i-1);
    end
end

for i = 1: length(PointCloud(:,1))
    PointCloud(i,:) = register(PointCloud(i,1),PointCloud(i,2),PointCloud(i,3),ImageOrientation,ImagePosition,PixelSpacing);
end
for i = 1: length(LandMarkS(:,1))
    LandMarkS(i,:) = register(LandMarkS(i,1),LandMarkS(i,2),LandMarkS(i,3),ImageOrientation,ImagePosition,PixelSpacing);
    LandMarkE(i,:) = register(LandMarkE(i,1),LandMarkE(i,2),LandMarkE(i,3),ImageOrientation,ImagePosition,PixelSpacing);
end

PointCloud(:,3) = - PointCloud(:,3);
LandMarkS(:,3)  = - LandMarkS(:,3);
LandMarkE(:,3)  = - LandMarkE(:,3);
%PointCloud(:,1) = - PointCloud(:,1);
%LandMarkS(:,1)  = - LandMarkS(:,1);
%LandMarkE(:,1)  = - LandMarkE(:,1);
%PointCloud(:,1) = - PointCloud(:,1);

%
figure(546)
scatter3(PointCloud(:,1),PointCloud(:,2),PointCloud(:,3));
hold on;
scatter3(LandMarkS(:,1),LandMarkS(:,2),LandMarkS(:,3),50,'r','filled');
scatter3(LandMarkE(:,1),LandMarkE(:,2),LandMarkE(:,3),50,'r','filled');
axis equal;

%XMeanApex = mean(PointCloud(1:SliceNum(1),1));
%YMeanApex = mean(PointCloud(1:SliceNum(1),2));

%% load stl from well-defined data
LAX2 = load(strcat(outputPath,caseName,'/ED','/LAX2-o.mat'),'-ascii');
LAX3 = load(strcat(outputPath,caseName,'/ED','/LAX3-o.mat'),'-ascii');
LAX4 = load(strcat(outputPath,caseName,'/ED','/LAX4-o.mat'),'-ascii');

for i = 1: length(LAX2(:,1))
    LAX2(i,:) = register(LAX2(i,1),LAX2(i,2),LAX2(i,3),ImageOrientation,ImagePosition,PixelSpacing);
end
for i = 1: length(LAX3(:,1))
   LAX3(i,:) = register(LAX3(i,1),LAX3(i,2),LAX3(i,3),ImageOrientation,ImagePosition,PixelSpacing);
end
for i = 1: length(LAX4(:,1))
    LAX4(i,:) = register(LAX4(i,1),LAX4(i,2),LAX4(i,3),ImageOrientation,ImagePosition,PixelSpacing);
end

LAX2(:,3) = - LAX2(:,3);
LAX3(:,3) = - LAX3(:,3);
LAX4(:,3) = - LAX4(:,3);
LAX2(:,1) = - LAX2(:,1);
LAX3(:,1) = - LAX3(:,1);
LAX4(:,1) = - LAX4(:,1);
%totalLAX = [LAX3;LAX4];
totalLAX = [LAX2;LAX3;LAX4];
LAXNum = length(totalLAX(:,1));
totalLAX = moveA2B(totalLAX,PointCloud);

%MinDisNum     =   [0.9000 -0.9000 0.9000];      % HYP01
%MinDelta      =   45;                           % HYP01
%MinDisNum     =   [0 0 0];                      % SC-HF-NI-33
%MinDelta      =   30;                           % SC-HF-NI-33

direction = [0 0 1];
centerX = mean(totalLAX(:,1));
centerY = mean(totalLAX(:,2));
origin = [centerX centerY 0];

%ToAdd = repmat(MinDisNum, LAXNum,1);
%totalLAX = totalLAX + ToAdd;
%totalLAX = rotateHori(totalLAX, direction, origin, MinDelta);

%LAX3 = totalLAX(1:length(LAX3(:,1)),:);
%LAX4 = totalLAX((length(LAX3(:,1))+1:(length(totalLAX))),:);

LAX2 = totalLAX(1:length(LAX2(:,1)),:);
LAX3 = totalLAX((length(LAX2(:,1))+1):(length(LAX2(:,1))+length(LAX3(:,1))),:);
LAX4 = totalLAX((length(LAX2(:,1))+length(LAX3(:,1))+1:(length(totalLAX))),:);


% read directly from the mat file
%LAX2 = load(strcat(outputPath,caseName,'/LAX2.mat'),'LAX2','-ascii');

%[normalLAX2, faceLAX2] = faceNormal(LAX2,1,20,40);
%%

%LAX3(:,3) = LAX3(:,3) -2 ;
%LAX4(:,3) = LAX4(:,3) -2 ;
%PointCloud(:,3) = PointCloud(:,3) -5;
sliceToAdj = 6;
offsetSize = -5;

%LAX2 = CorrectOffset(LAX2,faceLAX2,offsetSize);
%LAX2(:,2) = LAX2(:,2) - 10;
%LAX3(:,3) = LAX3(:,3)+10;
LAX3(:,1) = LAX3(:,1)-10;
LAX4(:,1) = LAX4(:,1)+5;
LAX4(:,2) = LAX4(:,2)+14;
LAX3(:,2) = LAX3(:,2)+5;
%LAX2(:,3) = LAX2(:,3)+5;
%PointCloud(:,3) = PointCloud(:,3) + 4;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
if(sliceToAdj == 1)
    PointCloud(1:SliceNumAll(sliceToAdj),:) = ...
    CorrectOffset(PointCloud(1:SliceNumAll(sliceToAdj),:),faceLAX2,offsetSize);
else
PointCloud(SliceNumAll(sliceToAdj-1)+1:SliceNumAll(sliceToAdj),:) = ...
    CorrectOffset(PointCloud(SliceNumAll(sliceToAdj-1)+1:SliceNumAll(sliceToAdj),:),faceLAX2,offsetSize);
end

% save(strcat(outputPath,caseName,'/Epi_0425.mat') , 'PointCloud', '-ascii');

%% plot here

h = figure(200)
%scatter3(PointCloud_i(:,1),PointCloud_i(:,2),PointCloud_i(:,3));
scatter3(PointCloud(:,1),PointCloud(:,2),PointCloud(:,3));
hold on;
xlabel('X');
scatter3(LandMarkS(:,1),LandMarkS(:,2),LandMarkS(:,3),50,'r','filled');
scatter3(LandMarkE(:,1),LandMarkE(:,2),LandMarkE(:,3),50,'r','filled');
scatter3(LAX4(:,1), LAX4(:,2), LAX4(:,3),'filled','r');     % yes
scatter3(LAX3(:,1), LAX3(:,2), LAX3(:,3),'filled','y');     % fan
scatter3(LAX2(:,1), LAX2(:,2), LAX2(:,3),'filled','g');     % fan
axis equal;
hold off;

%LAX2(:,3) = LAX2(:,3) + 4;
%LAX3(:,3) = LAX3(:,3) + 4;
%% plot 加密的 result
fMinResolution = 1;
PointCloudAll = [PointCloud;LAX2;LAX4];
%PointCloudAll = [PointCloud;LAX2;LAX3;LAX4];
SliceNumAll(length(SliceNumAll)+1) = length(LAX2(:,1)) + SliceNumAll(length(SliceNumAll));
%SliceNumAll(length(SliceNumAll)+1) = length(LAX3(:,1)) + SliceNumAll(length(SliceNumAll));
SliceNumAll(length(SliceNumAll)+1) = length(LAX4(:,1)) + SliceNumAll(length(SliceNumAll));

origin.vertices = PointCloudAll;
origin.faces = convhull(PointCloudAll(:,1),PointCloudAll(:,2),PointCloudAll(:,3));

[mesh.vertices, mesh.faces] = ...
               LoopSubdivisionLimited( origin.vertices, origin.faces, 2);
           
%[mesh.vertices, mesh.faces] = ...
%               LoopSubdivisionLimited( mesh.vertices, mesh.faces, fMinResolution);
%figure(246)
%patch(mesh,'FaceColor',[0.6 0.6 0.6],'facecolor','cyan','FaceAlpha',0.2,'EdgeAlpha',0.3);
           
%mesh.vertices=lpflow_trismooth(mesh.vertices,mesh.faces); 

% %计算法向量
[normal,~] = compute_normal(mesh.vertices, mesh.faces);
mesh.normal = normal';

% 修改法向量，使得指向mesh外部
%pmesh.faces = mesh.faces;
%pmesh.vertices = mesh.vertices;
%for i = 1:length(mesh.vertices(:,1))
%    i
%    endpoint = mesh.vertices(i,:) + mesh.normal(i,:).*20;
%  if(inpolyhedron(origin,endpoint))
%        mesh.normal(i,:) = -mesh.normal(i,:);
%    end
%end

%tmp.faces = mesh.faces;
%tmp.vertices = mesh.vertices;
%vnum = size(tmp.vertices(:,1));
%figure(233)
%   patch(tmp,'FaceColor',[0.6 0.6 0.6],'facecolor','cyan','FaceAlpha',0.2,'EdgeAlpha',0.3);
%   hold on;
%   quiver3(mesh.vertices(:,1),mesh.vertices(:,2),mesh.vertices(:,3),...
%       mesh.normal(:,1),mesh.normal(:,2),mesh.normal(:,3),'r');

    
    
    
%% attractor point
figure(233)
attractorPerAxis = 7;
deformSize = 0.1;
iteration = 5;
SkipSlice = [0];
pmesh = mesh;
 
%figure(456) 
for i = 1:5
    i
    %deformSize = -1.7 / i;
   attractor = getAttractor(SkipSlice,PointCloudAll,SliceNumAll,attractorPerAxis);
    [pmesh.deformScale, pmesh.diffVector] = calculateDeformFactor(pmesh,attractor);
    pmesh = deformDiscrepancy(pmesh, deformSize);
    
    %pmesh.vertices=lpflow_trismooth(pmesh.vertices,pmesh.faces); 
    
    patchmesh.vertices = pmesh.vertices;
    patchmesh.faces = pmesh.faces;
    clf;
    hold on;
    axis equal;
    patch(patchmesh,'FaceColor',[0.6 0.6 0.6],'facecolor','cyan','FaceAlpha',0.2,'EdgeAlpha',0.3);
    scatter3(attractor(:,1), attractor(:,2), attractor(:,3),'filled','b'); % yes
    %quiver3(patchmesh.vertices(:,1),patchmesh.vertices(:,2),patchmesh.vertices(:,3),...
    %    pmesh.diffVector(:,1),pmesh.diffVector(:,2),pmesh.diffVector(:,3),'r');
    %   quiver3(patchmesh.vertices(1:2:vnum,1),patchmesh.vertices(1:2:vnum,2),patchmesh.vertices(1:2:vnum,3),...
 %       pmesh.diffVector(1:2:vnum,1),pmesh.diffVector(1:2:vnum,2),pmesh.diffVector(1:2:vnum,3),'r');
%    scatter3(PointCloud(:,1),PointCloud(:,2),PointCloud(:,3));
%    scatter3(LAX2(:,1), LAX2(:,2), LAX2(:,3),'y'); % yes
%    scatter3(LAX3(:,1), LAX3(:,2), LAX3(:,3),'g'); % yes
%    scatter3(LAX4(:,1), LAX4(:,2), LAX4(:,3),'r'); % yes
    lighting gouraud;set(gcf,'color','white');
    %title(strcat('aligned ToPoint deformSize = 0.06; iteration = 30; MeanDiscrepance = ',num2str(mean(abs(distances)))));
    %view(45,45);
    view(90,0);
    pause(0.1);
end

pmesh.vertices=lpflow_trismooth(pmesh.vertices,pmesh.faces); 
[pmesh.vertices, pmesh.faces] = ...
               LoopSubdivisionLimited( pmesh.vertices, pmesh.faces, 2);

           
           
%% 把中轴掰直的过程

tmp_mid = pmesh.vertices;
p1 = LandMarkS(2,:); %Apex
p2 = LandMarkE(2,:); %Basal
% find upper and lower point here
ind = find(pmesh.vertices(:,3) == min(pmesh.vertices(:,3)));
Basal = pmesh.vertices(ind,:);


%Apex = [183, 173, 72];
ROT = p1 - p2;
X = ROT(1);
Y = ROT(2);
Z = ROT(3);

%Psi = acos(Y/(sqrt(X^2 + Y^2)))*180/pi;
Psi = acos(Y/(sqrt(X^2 + Y^2)))*180/pi;

%Psi = 180;
RzNPsi = [  cosd(-Psi)  -sind(-Psi) 0;
            sind(-Psi)  cosd(-Psi)  0;
            0           0           1];

pointNum = length(tmp_mid(:,1));
origin = repmat(p2, pointNum,1);

figure(233)
hold on;
axis equal;
set(gca,'XGrid','on');
%axis([-100 200 -100 200 -100 100])
quiver3(p2(1),p2(2),p2(3),ROT(:,1),ROT(:,2),ROT(:,3),'r');
scatter3(tmp_mid(1:10:end,1),tmp_mid(1:10:end,2),tmp_mid(1:10:end,3));
scatter3(Basal(:,1),Basal(:,2),Basal(:,3),'b','filled');
scatter3(p1(:,1),p1(:,2),p1(:,3),'b','filled');
scatter3(p2(:,1),p2(:,2),p2(:,3),'b','filled');
xlabel('X');
ylabel('Y');
zlabel('Z');

Basal = Basal - origin(1,:);
Basal = Basal * RzNPsi;
Basal = Basal + origin(1,:);
originBasal = repmat(Basal, pointNum,1);

tmp_mid = (tmp_mid - origin) * RzNPsi + origin - originBasal;

Basal = Basal - originBasal(1,:);

p1 = p1 - origin(1,:);
p1 = p1 * RzNPsi;
p1 = p1 + origin(1,:) - originBasal(1,:);

p2 = p2 - origin(1,:);
p2 = p2 * RzNPsi;
p2 = p2 + origin(1,:) - originBasal(1,:);

%tmp_mid(:,1) = - tmp_mid(:,1);
%Basal(:,1) = - Basal(:,1);
%p1(:,1) = - p1(:,1);
%p2(:,1) = - p2(:,1);

ROT = p1 - p2;
figure(22)
    hold on;
    axis equal;
    set(gca,'XGrid','on');
    axis([])
    quiver3(p2(1),p2(2),p2(3),ROT(:,1),ROT(:,2),ROT(:,3),'r');
    %patch(deleteme,'FaceColor',[0.6 0.6 0.6],'facecolor','cyan','FaceAlpha',0.2,'EdgeAlpha',0.3);
    scatter3(tmp_mid(1:10:end,1),tmp_mid(1:10:end,2),tmp_mid(1:10:end,3));
    scatter3(Basal(:,1),Basal(:,2),Basal(:,3),'g','filled');
    scatter3(p1(:,1),p1(:,2),p1(:,3),'b','filled');
    scatter3(p2(:,1),p2(:,2),p2(:,3),'b','filled');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    hold off;

%% save result
pmesh.vertices = tmp_mid;
pmesh.p1 = p1;
pmesh.p2 = p2;
pmesh.Basal = Basal;

%
outputP = strcat('./output/0518_landmark_aligned/EpiED/mesh/new0528/',caseName);
mkdir(strcat(outputP));
save(strcat(outputP,'/pmesh.mat'),'pmesh');
%save(strcat(outputP,'/LAX3.mat'),'LAX3');
%save(strcat(outputP,'/LAX4.mat'),'LAX4');
save(strcat(outputP,'/SliceNumAll.mat'),'SliceNumAll','-ascii');
save(strcat(outputP,'/SliceNum.mat'),'SliceNum','-ascii');

%% plot the result
tmp.vertices = pmesh.vertices;
tmp.faces = pmesh.faces;
figure(240)  
hold on;
axis equal;
patch(tmp,'FaceColor',[0.6 0.6 0.6],'facecolor','cyan','FaceAlpha',0.2,'EdgeAlpha',0.3);
%scatter3(PointCloud(:,1),PointCloud(:,2),PointCloud(:,3));
scatter3(Basal(1),Basal(2),Basal(3),100,'r','filled');





