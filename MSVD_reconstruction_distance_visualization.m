clear all;
caseList = {{
'SC-HF-I-01'
'SC-HF-I-02'
'SC-HF-I-04'
'SC-HF-I-05'
'SC-HF-I-06'
'SC-HF-I-11'
'SC-HF-I-12'
'SC-HF-I-40'}
{
'SC-HF-NI-03'
'SC-HF-NI-04'
'SC-HF-NI-07'
'SC-HF-NI-12'
'SC-HF-NI-13'
'SC-HF-NI-14'
'SC-HF-NI-33'
'SC-HF-NI-34'}
{
'SC-HYP-01'
'SC-HYP-03'
'SC-HYP-06'
'SC-HYP-07'
'SC-HYP-08'
'SC-HYP-09'
'SC-HYP-10'
'SC-HYP-11'}
};

rootpath = 'C:\Users\37908\Desktop\tmp\mesh_recons_0921\normalized\Epi\shell_nii\';
deformed_root= 'C:\Users\37908\Desktop\tmp\normalized\result_17_aod_922\';
listt = dir(deformed_root);

% parameter
highlight_percentage = 0.02;
testQ = 80; maxK = 5;
% parameter

caseType = 2;
TXX = zeros(50,50,50,8);
for tt = 1:8
    tt
    curName = caseList{caseType}{tt};
    F = load([deformed_root, curName,'\all_vari']);
    mesh = F.epi; 
    mesh = unique(round(mesh),'rows');
    TXX(:,:,:,tt) = pc_2_nii(mesh);
end

%%
for caseIdx = 4:8
    close;
    label = load(['C:\Users\37908\Desktop\zzz\', 'label_3_24.txt']);
    len = sum(label==caseType);
    gndTX = linspace(caseType,caseType,len);

    meanTX = mean(TXX,4);
    meanTXM = repmat(meanTX,1,1,1,len);
    TX = TXX - meanTXM;

    % MPCA
     disp('----------Calculating MPCA----------')
     [tUs, odrIdx, ~, ~]  = MPCA(TXX,gndTX,testQ,maxK);
     disp('----------------Done----------------')

    %%
    crtCaseAp = TXX(:,:,:,caseIdx);
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

    res = abs(res) + meanTX;
    res_p = abs(res);
    top_res_p = top_res;

    % 负值变0
    res_p(res<0) = 0;
    top_res_p(top_res<0) = 0;

    % 线性变换到[0 255]
    res_plot = (res_p-min(min(min(res_p))))/(max(max(max(res_p))))*255;
    top_res_plot = (top_res-min(min(min(top_res_p))))/(max(max(max(res_p))))*255;

    %% With top features shown :3D Visualization

    threshhold = 5;
    all_plot = to_plot_point_cloud(res_plot, threshhold);
    top_plot = to_plot_point_cloud(top_res_plot, threshhold);

    figure(234)
        subplot(121)
            scatter3(all_plot(:,1),all_plot(:,2),all_plot(:,3),all_plot(:,4)/2,'MarkerEdgeColor','none','MarkerFaceColor',[0 .75 .75]);
            hold on;
            scatter3(top_plot(:,1),top_plot(:,2),top_plot(:,3),top_plot(:,4)*5,'MarkerEdgeColor','none','MarkerFaceColor',[1 0 0]);
            axis equal
            xlabel('X axis');
            ylabel('Y axis');
            zlabel('Z axis');
            axis([min(all_plot(:,1)) max(all_plot(:,1)) min(all_plot(:,2)) max(all_plot(:,2)) min(all_plot(:,3)) max(all_plot(:,3))])
            title([caseName,'-Epi-Reconstructed']);
pause()
  %      subplot(122)
  %          patch(origin,'FaceColor',[0.6 0.6 0.6],'facecolor','cyan','FaceAlpha',0.2,'EdgeAlpha',0.15);
  %          axis([-40 40 -40 40 0 100])
  %          axis equal
  %          xlabel('X axis');
  %          ylabel('Y axis');
  %          zlabel('Z axis');
  %          title([caseName,'-Epi-Origin']);
    
end