%% 
% input: all dicom and txt files of this case
% output: the mesh construction result of this case, including vertices, faces,discrepancy, confidence interval  

% process: 
%    read the xyz coordinate of contours from the txt file
%    transform the xyz coordinate to rlapfh depending on the information of corresponding dcm file
%    add the LAX contour to the vertices
%    register the pointcloud
%    subdivision
%    use attractor to perform deformation
%    make the Z-axis of the mesh vertical
%    output the mesh information in the form of mat file 
    
    
%
%%
clear all;
caseName = 'SC-HF-I-06';

txtPath = strcat('./data/',caseName,'/txt/contours-manual/IRCCI-expert/');
dcmPath = strcat('./data/',caseName,'/dcm/'); 
LAXFilePath = strcat('./data/',caseName,'/LAX/'); 
outputMatPath = strcat('./data/',caseName,'/output/'); 

% choose to construct the mesh of Endo/Epi ES/ED
readEndoES = false;
readEndoED = false;
readEpiES = true;

%%
% read txt file
flagi = 0;flago = 0;

listing = dir(txtPath);
numofFile = length(listing);
PointCloud = [];
SliceNum = [];
SliceNumTotal = [];

dcmList = dir(dcmPath);
Corresponding = dcmList(3).name(1:7); 

numO = 0;
for i = 3:numofFile
% for i = 5:5
    fileName = listing(i).name;
    ct = load(strcat(txtPath,fileName));
    IDString = fileName(9:12);
    ID = str2double(fileName(10:12));
    numSlice = ceil(ID/20);
    fileType = fileName(14);
    % if the points represent endocardial 
    if(fileType == 'i' && ~(flagi == numSlice) && readEndoES)
        numO = numO + 1;
        SliceNum(numO) = length(ct(:,1));
        type = 1;newSliceGet = true;flagi = numSlice;
        c = linspace(1,1,length(ct(:,1)));
        xi = ct(:,1)';
        yi = ct(:,2)';
        zi = c';
        fprintf('Endo - ES - ID is %s\n',IDString);
    elseif(fileType == 'i' && (flagi == numSlice) && readEndoED)
        numO = numO + 1;
        SliceNum(numO) = length(ct(:,1));
        type = 2;newSliceGet = true;
        c = linspace(1,1,length(ct(:,1)));
        xied = ct(:,1)';
        yied = ct(:,2)';
        zied = c';
        fprintf('Endo - ED - ID is %s\n',IDString);
    % else if the points represent epicardium  
    elseif(fileType == 'o' && ~(flago == numSlice) && readEpiES)
        numO = numO + 1;
        SliceNum(numO) = length(ct(:,1));
        type = 3;newSliceGet = true;flago = numSlice;
        c = linspace(1,1,length(ct(:,1)));
        xo = ct(:,1)';
        yo = ct(:,2)';
        zo = c';
        fprintf('Epi - ED - ID is %s\n',IDString);
    else
        newSliceGet = false;
    end
    
    % read corresponding dcm
    if(newSliceGet)
        %EndoESPoints = [xi',yi',zi'];
        %EpiESPoints = [xo',yo',zo'];
        dcmName = strcat(Corresponding,'-',IDString,'.dcm');
        DCM = dicominfo(strcat(dcmPath,dcmName));
        % get information from the dcm file
        ImageOrientation        =    DCM.ImageOrientationPatient';
        PixelSpacing            =    DCM.PixelSpacing;
        SliceThickness          =    DCM.SliceThickness;
        ImagePosition           =    DCM.ImagePositionPatient';
        
        % transform this txt file, keep the position
        if(type == 1)
            for idx = 1 : length(xi)
                newPos = xyz2rlapfh(xi(idx),yi(idx),zi(idx),ImageOrientation,ImagePosition,PixelSpacing,SliceThickness);
                PointCloud = [PointCloud; newPos];
            end
            
        elseif(type == 3)
            % epi es
            for idx = 1 : length(xo)
                newPos = xyz2rlapfh(xo(idx),yo(idx),zo(idx),ImageOrientation,ImagePosition,PixelSpacing,SliceThickness);
                PointCloud = [PointCloud; newPos];
            end
        elseif(type == 2)
            % endo ed
            for idx = 1 : length(xied)
                newPos = xyz2rlapfh(xied(idx),yied(idx),zied(idx),ImageOrientation,ImagePosition,PixelSpacing,SliceThickness);
                PointCloud = [PointCloud; newPos];
            end
        end
    end
end


%figure(546)
%scatter3(PointCloud(:,1),PointCloud(:,2),PointCloud(:,3));
%% register

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
PointCloud(:,3) = - PointCloud(:,3);

XMeanApex = mean(PointCloud(1:SliceNum(1),1));
YMeanApex = mean(PointCloud(1:SliceNum(1),2));

%% load LAX information stl from well-defined data

LAX2 = load(strcat(LAXFilePath,caseName,'/LAX2-o.mat'),'LAX2','-ascii');
LAX3 = load(strcat(LAXFilePath,caseName,'/LAX3-o.mat'),'LAX3','-ascii');
LAX4 = load(strcat(LAXFilePath,caseName,'/LAX4-o.mat'),'LAX4','-ascii');

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

totalLAX = [LAX2;LAX3;LAX4];
LAXNum = length(totalLAX(:,1));
totalLAX = moveA2B(totalLAX,PointCloud);

LAX2 = totalLAX(1:length(LAX2(:,1)),:);
LAX3 = totalLAX((length(LAX2(:,1))+1):(length(LAX2(:,1))+length(LAX3(:,1))),:);
LAX4 = totalLAX((length(LAX2(:,1))+length(LAX3(:,1))+1:(length(totalLAX))),:);

%% plot here

h = figure(233)
%scatter3(PointCloud_i(:,1),PointCloud_i(:,2),PointCloud_i(:,3));
scatter3(PointCloud(:,1),PointCloud(:,2),PointCloud(:,3));
hold on;axis equal;
scatter3(LAX2(:,1), LAX2(:,2), LAX2(:,3),'filled','y');  % fan
scatter3(LAX3(:,1), LAX3(:,2), LAX3(:,3),'filled','g'); % yes
scatter3(LAX4(:,1), LAX4(:,2), LAX4(:,3),'filled','r'); % yes
hold off;


%% use subdivision to make the vertices denser
fMinResolution = 1;
PointCloudAll = [PointCloud;LAX2;LAX3;LAX4];
SliceNumAll(length(SliceNumAll)+1) = length(LAX2(:,1)) + SliceNumAll(length(SliceNumAll));
SliceNumAll(length(SliceNumAll)+1) = length(LAX3(:,1)) + SliceNumAll(length(SliceNumAll));
SliceNumAll(length(SliceNumAll)+1) = length(LAX4(:,1)) + SliceNumAll(length(SliceNumAll));

origin.vertices = PointCloudAll;
origin.faces = convhull(PointCloudAll(:,1),PointCloudAll(:,2),PointCloudAll(:,3));

[mesh.vertices, mesh.faces] = ...
               LoopSubdivisionLimited( origin.vertices, origin.faces, 2);
           
%% attractor point
figure(233)
attractorPerAxis = 7;
deformSize = 0.06;
IterationTime = 10;
SkipSlice = [0];
pmesh = mesh;
 
% perform attractor for 10 times 
for i = 1:IterationTime
    fprintf('%d\n',i);
    %deformSize = 0.3 / i;
    attractor = getAttractor(SkipSlice,PointCloudAll,SliceNumAll,attractorPerAxis);
    [pmesh.deformScale, pmesh.diffVector] = calculateDeformFactor(pmesh,attractor);
    pmesh = deformDiscrepancy(pmesh, deformSize); 
    
    patchmesh.vertices = pmesh.vertices;
    patchmesh.faces = pmesh.faces;
    clf;
    hold on;
    axis equal;
    patch(patchmesh,'FaceColor',[0.6 0.6 0.6],'facecolor','cyan','FaceAlpha',0.2,'EdgeAlpha',0.3);
    scatter3(attractor(:,1), attractor(:,2), attractor(:,3),'filled','b'); % yes
    view(90,0);
    pause(1);
end

% subdivide again
pmesh.vertices=lpflow_trismooth(pmesh.vertices,pmesh.faces); 
[pmesh.vertices, pmesh.faces] = ...
               LoopSubdivisionLimited( pmesh.vertices, pmesh.faces, 2);

           
%% 把中轴掰直的过程
tmp = pmesh.vertices;

% find upper and lower point here
ind = find(pmesh.vertices(:,3) == min(pmesh.vertices(:,3)));
Basal = pmesh.vertices(ind,:);
%XMeanApex
ind = find(pmesh.vertices(:,3) == max(pmesh.vertices(:,3)));
TopZ = pmesh.vertices(ind,3);
Apex = [XMeanApex, YMeanApex, TopZ];

ROT = Apex - Basal;
X = ROT(1);
Y = ROT(2);
Z = ROT(3);

Psi = acos(Y/(sqrt(X^2 + Y^2)))*180/pi;
Phi = acos(Z/(sqrt(X^2 + Y^2 + Z^2)))*180/pi;

x = sind(Phi)*cosd(Psi);
y = sind(Phi)*sind(Psi);
z = cosd(Phi);

RzNPsi = [  cosd(-Psi)  -sind(-Psi) 0;
            sind(-Psi)  cosd(-Psi)  0;
            0           0           1];

RxNPxi = [   1      0           0 ;
             0      cosd(-Phi)  -sind(-Phi);
             0      sind(-Phi)  cosd(-Phi)];

pointNum = length(tmp(:,1));
origin = repmat(Basal, pointNum,1);
 
figure(234)
hold on;
axis equal;
set(gca,'XGrid','on');
axis([100 200 100 200 -100 100])
quiver3(Basal(1),Basal(2),Basal(3),ROT(:,1),ROT(:,2),ROT(:,3),'r');
scatter3(tmp(1:10:end,1),tmp(1:10:end,2),tmp(1:10:end,3));
scatter3(Basal(:,1),Basal(:,2),Basal(:,3),'b','filled');
scatter3(Apex(:,1),Apex(:,2),Apex(:,3),'b','filled');

tmp = tmp - origin;
tmp = tmp * RzNPsi * RxNPxi;

Apex_Ali = Apex - origin(1,:);
Apex_Ali = Apex_Ali * RzNPsi * RxNPxi;

Basal_Ali = Basal - origin(1,:);
Basal_Ali = Basal_Ali * RzNPsi * RxNPxi;

figure(235)
    hold on;
    axis equal;
    set(gca,'XGrid','on');
    %axis([])
    %quiver3(Basal(1),Basal(2),Basal(3),ROT(:,1),ROT(:,2),ROT(:,3),'r');
    %scatter3(tmp(:,1),tmp(:,2),tmp(:,3));
    scatter3(tmp(1:10:end,1),tmp(1:10:end,2),tmp(1:10:end,3));
    scatter3(Basal_Ali(:,1),Basal_Ali(:,2),Basal_Ali(:,3),'g','filled');
    scatter3(Apex_Ali(:,1),Apex_Ali(:,2),Apex_Ali(:,3),'g','filled');
    hold off;
    
pmesh.vertices = tmp;

%% calculate the discrepancy, confidence interval and then plot

toPlot.faces =pmesh.faces;
toPlot.vertices = pmesh.vertices;
vnum = size(toPlot.vertices(:,1));

[distances,~,~,~,~] = point2trimesh(toPlot, 'QueryPoints', attractor); 
[~,~,muci_before] = normfit(abs(distances'),0.25);

discrepancy = mean(abs(distances));

figure(221) 
    clf;
    hold on;
    axis equal;
    patch(toPlot,'FaceColor',[0.6 0.6 0.6],'facecolor','cyan','FaceAlpha',0.2,'EdgeAlpha',0.3);
    Tit1 = strcat('DeformSize = ',num2str(deformSize),'; Iteration = ',num2str(IterationTime),'; MeanDiscrepance =',num2str(discrepancy));
    Tit2 = strcat('Vertice num = ',num2str(vnum(1)),';75% Confidence interval = [',num2str(muci_before(1)),', ',num2str(muci_before(2)),']');
    %Tit = strcat('DeformSize = 0.06; Iteration = 30; MeanDiscrepance =',num2str(discrepancy),';vertice num = ',num2str(vnum(1)),';75% Confidence interval = [',num2str(muci_before(1)),', ',num2str(muci_before(2)),']');
    title({Tit1;Tit2})
    view(45,45);
    
%% save result
mkdir(strcat(outputMatPath,caseName));
save(strcat(outputMatPath,caseName,'/pmesh.mat'),'pmesh');
