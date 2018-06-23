function surface_points_before = Input_TXTName_Output_PointsAfterButterflySubdivision(fileName)


    fileName = 'SC-HF-I-01';
    LaplacianSmoothTimes = 1;
    attractorPointDensity = 4;
    % 8 is read from the dcm file of this case
    SliceHeight = 8;

    patientPathIn = '/contours-manual/IRCCI-expert/';
    dataPath = './data/dataset/';

    outputFile = '0409Contours_Subdivision/';

    listingP = dir(strcat(dataPath,fileName,patientPathIn));
    numofPatient = length(listingP);

    xi = [];
    yi = [];
    zi = [];
    xo = [];
    yo = [];
    zo = [];
    
    flagi = 0;
    flago = 0;
    flagp = 0;
    
    numofFile = length(listingP);
    
    % for i = 3:3
    for i = 3:numofFile
        contourName = listingP(i).name;
        % if(fileName == )
        ct = load(strcat(dataPath, fileName, patientPathIn, contourName));
        height = str2double(contourName(10:12));
        numSlice = ceil(height/20);
        if(rem(height, 20)==0)
            numFrame = 20;
        else
            numFrame = rem(height, 20);
        end
        heightSlice = ceil(height/20);
        heightArr = linspace(heightSlice*SliceHeight,heightSlice*SliceHeight,length(ct(:,1)));
        frameArr = linspace(numFrame,numFrame,length(ct(:,1)));
        
        fileType = contourName(14);
        % if the points represent endocardial 
        if(fileType == 'i' && ~(flagi == numSlice))
            flagi = numSlice;
            c = linspace(1,1,length(ct(:,1)));
            xi = [xi ct(:,1)'];
            yi = [yi ct(:,2)'];
            zi = [zi heightArr];
        %    [xiAft,yiAft,ziAft] = Subdivision_Points_Coor_0409([xi',yi',zi']);
        % else if the points represent epicardium  
        elseif(fileType == 'o' && ~(flago == numSlice))
            flago = numSlice;
            c = linspace(120,120,length(ct(:,1)));
            xo = [xo ct(:,1)'];
            yo = [yo ct(:,2)'];
            zo = [zo heightArr];
        end
    end

    % the 3d information about endocardial and epicardium
    EndocardiumBeforeSmooth = [xi',yi',zi'];
    EpicardiumBeforeSmooth = [xo',yo',zo'];
  
    % add upper and lower limit point the the point cloud
    EndocardiumBeforeSmoothAdded = addUpandDownLimit(EndocardiumBeforeSmooth);
    EpicardiumBeforeSmoothAdded = addUpandDownLimit(EpicardiumBeforeSmooth);

    % Laplacian smooth
    EndocardiumAfterSmooth = LaplacianSmooth(EndocardiumBeforeSmooth, LaplacianSmoothTimes);
    EpicardiumAfterSmooth = LaplacianSmooth(EpicardiumBeforeSmooth, LaplacianSmoothTimes);

    % add upper and lower limit point the the point cloud
    EndocardiumAfterSmoothAdded = addUpandDownLimit(EndocardiumAfterSmooth);
    EpicardiumAfterSmoothAdded = addUpandDownLimit(EpicardiumAfterSmooth);
    
    
%% analysis
    x_before = EpicardiumBeforeSmoothAdded(:,1);
    y_before = EpicardiumBeforeSmoothAdded(:,2);
    z_before = EpicardiumBeforeSmoothAdded(:,3);
    faces_before = convhull(x_before,y_before,z_before,'simplify',true);

    x_after = EndocardiumAfterSmoothAdded(:,1);
    y_after = EndocardiumAfterSmoothAdded(:,2);
    z_after = EndocardiumAfterSmoothAdded(:,3);
    faces_after = convhull(x_after,y_after,z_after,'simplify',true);

   
    %% Bufferfly subdivision
        
        fMinResolution = 0.01;
        [mfRefinedMeshBefore, mnTriangulationBefore] = ...
           LoopSubdivisionLimited( EpicardiumBeforeSmoothAdded, faces_before, fMinResolution);
        x_before = mfRefinedMeshBefore(:,1);
        y_before = mfRefinedMeshBefore(:,2);
        z_before = mfRefinedMeshBefore(:,3);


        [mfRefinedMeshAfter, mnTriangulationAfter] = ...
           LoopSubdivisionLimited( EndocardiumAfterSmoothAdded, faces_after, fMinResolution);
        x_after = mfRefinedMeshAfter(:,1);
        y_after = mfRefinedMeshAfter(:,2);
        z_after = mfRefinedMeshAfter(:,3);


        %% farthest point sampling 

        FV_before.faces = mnTriangulationBefore;
        FV_before.vertices = mfRefinedMeshBefore;
        initialPoints = EpicardiumBeforeSmooth(1:4:length(EpicardiumBeforeSmooth(:,1)),:);

        FV_after.faces = mnTriangulationAfter;
        FV_after.vertices = mfRefinedMeshAfter;
        % initialPoints = EndocardiumBeforeSmooth(1:1:length(EndocardiumBeforeSmooth(:,1)),:);

        [~,surface_points_before] = point2trimesh(FV_before, 'QueryPoints', initialPoints); 
        [~,surface_points_after] = point2trimesh(FV_after, 'QueryPoints', initialPoints); 
        
        
       % make the height of surface_points_before integer
       surface_points_before(:,3) = initialPoints(:,3);
