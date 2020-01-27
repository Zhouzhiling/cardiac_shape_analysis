clear
root = 'F:\Harvard\shape_analysis_before_1009\tmp\normalized_no\';
inroot = [root,'Input\'];
outroot = [root,'Output\'];
sourcepath = {[root,'Input\Endo\Template.nii'],[root,'Input\Epi\Template.nii']};

endosuffix = 'Endo\';
episuffix = 'Epi\';
suffix = {episuffix};%,endosuffix};
endofile = dir([inroot,endosuffix,'*']);
endonum = length(endofile);
epifile = dir([inroot,episuffix,'*']);
epinum = length(epifile);

ptnum=50;
TXX = zeros(ptnum,ptnum,ptnum,8);

[cover_height, cover_area] = calculate_cover_height_area(sourcepath);
chosenGroup = 'SC-HF-I';
count = 0;
for t = 3:27
    %t = 3;
    caseList = endofile(t).name;%;{'SC-HF-NI-03.nii'}%
    caseNameOnly = caseList(1:end-4);
%    if(~contains(caseNameOnly,chosenGroup))
%        continue;
%    end
    count = count + 1;
    disp(caseNameOnly);
	NAME{count} = caseNameOnly;
    mesh = cell(2,2);
    for i = 1:length(suffix)
        %close;
        caseName = caseList(1:end-4);
        cur_case = read_nii([inroot,suffix{i},caseName,'.nii']);
        cur_case_pc = to_plot_point_cloud(cur_case,10);
      %  targetpath = [root,'Input_0808_all\',caseName,'.nii'];

        vfxpath = [root,'Output\',suffix{i},caseName,'-_VelocityField_X.nii'];
        vfypath = [root,'Output\',suffix{i},caseName,'-_VelocityField_Y.nii'];
        vfzpath = [root,'Output\',suffix{i},caseName,'-_VelocityField_Z.nii'];
        aodpath = [root,'Output\',suffix{i},caseName,'-_TotalAOD.nii'];
        finaldefpath = [root,'Output\',suffix{i},caseName,'-_FinalDefSrc.nii'];
        
        source = read_nii(sourcepath{i});
        source_pc = to_plot_point_cloud(source,10);
      %  target = read_nii(targetpath);
        % velocity field x
        vfx_ = read_nii(vfxpath);
        vfy_ = read_nii(vfypath);
        vfz_ = read_nii(vfzpath);
        final_def = read_nii(finaldefpath);
        
        vfx = vfx_(:,:,:,end);
        vfy = vfy_(:,:,:,end);
        vfz = vfz_(:,:,:,end);
        
       % vfx = sum(vfx_,4);
       % vfy = sum(vfy_,4);
       % vfz = sum(vfz_,4);
        
        scale_aod = sqrt(vfx .^2 + vfy .^2 + vfz .^2);
        scale_aod(source==0)=0;
        TXX(:,:,:,count) = scale_aod;
%        tp = to_plot_point_cloud(scale_aod,0.2);
%        figure(1646)
%            scatter3(tp(:,1),tp(:,2),tp(:,3));
%            title(caseName);
%            pause
    end
end


%%
sz = size(TXX);
numSpl = sz(end);
N = length(sz)-1;

chosen_class = 9:16; 
feature = TXX(:,:,:,[chosen_class 25]);
mean_TXX = mean(feature,N+1);
fea3Dctr = feature-repmat(mean_TXX,[ones(1,N), 9]);%Centering
TX = fea3Dctr;
gnd = [linspace(0,0,8),1];

MPCADADim   =378;  % 200 most discriminative MPCA features for LDA
testQ       =90;   % Keep 97% variation in each mode
maxK        =5;    % One iteration only
[tUs,odrIdx,TXmean,Wgt,LDAU] = MPCALDA(feature,gnd,MPCADADim,testQ,maxK);

%%
outroot = 'F:\Document_Matlab\0113_for_ISBI\MPCA_HighLight\';
highlight_percentage = 0.08;
for caseIdx = [7]
%caseIdx = 3
    crtCaseAp = feature(:,:,:,caseIdx);
    Mat = TX(:,:,:,caseIdx);

    % TX是已经减过平均的
    xtus = tUs{1};
    ytus = tUs{2};
    ztus = tUs{3};

    length_x = size(xtus,1);
    length_y = size(ytus,1);
    length_z = size(ztus,1);

    featureNum = length_x * length_y * length_z;
    % for fidx = 1:featureNum
    res = zeros(50,50,50);
    top_res = zeros(50,50,50);
    for fidx = 1:featureNum
        disp(fidx);
        odr = odrIdx(fidx);

        idx_z = ceil(odr/(length_x * length_y));
        odr__ = odr - (idx_z-1)*(length_x * length_y);
        idx_y = ceil(odr__/length_x);
        idx_x = rem(odr__,length_x);
        if(idx_x==0) 
            idx_x = length_x ;
        end

        vector_x = xtus(idx_x,:)';
        vector_y = ytus(idx_y,:)';
        vector_z = ztus(idx_z,:)';
        sz_Mat = size(Mat);
        % 以下实现 res = vector_x' * vector_x * X + mean;
        tMat = zeros(sz_Mat);
        rrrx = zeros(sz_Mat);
        for t = 1:sz_Mat(3)
            tMat(:,:,t) = Mat(:,:,t)';
        end
        rx = mmat(tMat,vector_x);
        rrx = mmat(rx,vector_x');
        for t = 1:sz_Mat(3)
            rrrx(:,:,t) = rrx(:,:,t)';
        end
        res_x = rrrx;

        ry = mmat(Mat,vector_y);
        rry = mmat(ry,vector_y');
        res_y = rry;

        C_z = reshape(reshape(Mat,[],sz_Mat(3)) * vector_z,[sz_Mat(1:2) 1]);
        CC_z = reshape(reshape(C_z,[],1) * vector_z', sz_Mat);
        res_z = CC_z;
        res = res + res_x + res_y + res_z;
        if(fidx < highlight_percentage * featureNum)
            top_res = top_res + res_x + res_y + res_z;
        end
    end

    res = abs(res) + mean_TXX;
    res_p = abs(res);
    top_res_p = top_res;

    % 负值变0
    res_p(res<0) = 0;
    top_res_p(top_res<0) = 0;

    % 线性变换到[0 255]
    res_plot = (res_p-min(min(min(res_p))))/(max(max(max(res_p))))*255;
    top_res_plot = (top_res-min(min(min(top_res_p))))/(max(max(max(res_p))))*255;

    threshhold = 17;
    all_plot = to_plot_point_cloud(res_plot, threshhold);
    top_plot = to_plot_point_cloud(top_res_plot, 2);
    
    highlight_idx = find_cores_pt(all_plot,top_plot);
    
    % 黑白
%    clr = zeros(length(highlight_idx),3);
%    clr(highlight_idx==1,1:3) = 0.7;
    
    % 红色
    clr = zeros(length(highlight_idx),3);
    clr(:,2:3) = 0.5;
    clr(highlight_idx==1,2:3) = 0.8;
        
    % find another way to visualize MPCA result
    MESH.faces = convhull(all_plot(:,1:3));
    MESH.vertices = all_plot(:,1:3);
    
    MESH.vertices = lpflow_trismooth(MESH.vertices,MESH.faces); 
    MESH.vertices = lpflow_trismooth(MESH.vertices,MESH.faces); 
    MESH.vertices = lpflow_trismooth(MESH.vertices,MESH.faces); 
    
%     figure(233)
%         hold on;axis equal;xlabel('X');ylabel('Y'),zlabel('Z');
%         scatter3(all_plot(:,1),all_plot(:,2),all_plot(:,3),all_plot(:,4)/2,'MarkerEdgeColor','none','MarkerFaceColor',[0 .75 .75]);
%         scatter3(top_plot(:,1),top_plot(:,2),top_plot(:,3),top_plot(:,4)*10,'r','filled');
%         title(num2str(caseIdx));
%     
	close all;
    figure(234)
 %       camlight('headlight')
        hold on;axis equal;xlabel('X');ylabel('Y'),zlabel('Z');
  %       scatter3(all_plot(:,1),all_plot(:,2),all_plot(:,3),all_plot(:,4)/2,'MarkerEdgeColor','none','MarkerFaceColor',[0 .75 .75]);
        patch('Faces',MESH.faces,'Vertices',MESH.vertices,'FaceVertexCData',clr,'FaceColor','interp','EdgeAlpha',0);
		view([85,10]);
	outpath = [outroot,NAME{caseIdx+8},'.fig'];
	saveas(gcf,outpath);
	outpath2 = [outroot,NAME{caseIdx+8},'.png'];
	saveas(gcf,outpath2);
 %   [MESH,highlight_idx] = delete_useless_pt(MESH,highlight_idx);
%    looptime = 0;
%    
%    while(length(MESH.vertices)<3000 && looptime<5)
%        looptime = looptime + 1;
%        [MESH.vertices, MESH.faces] = ...
%                   LoopSubdivisionLimited(MESH.vertices, MESH.faces, 1);
end
        
    
    