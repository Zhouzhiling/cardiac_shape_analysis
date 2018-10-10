clear;
root = 'C:\Users\37908\Desktop\tmp\';
inroot = [root,'Input\'];
outroot = [root,'Output\'];
resroot = [root,'result_17_aod_921\'];
sourcepath = {[root,'Input\EndoES_Template.nii'],[root,'Input\EpiES_Template.nii']};
       
endosuffix = 'Endo\';
episuffix = 'Epi\';
suffix = {endosuffix,episuffix};
endofile = dir([inroot,endosuffix,'*']);
endonum = length(endofile);
epifile = dir([inroot,episuffix,'*']);
epinum = length(epifile);

for t = 3:epinum
    t
    %t = 3;
    caseList = {endofile(t).name;epifile(t).name};
    caseNameOnly = caseList{1}(1:end-4);
    respath = [resroot,caseNameOnly,'\'];
    mkdir(respath);
    mesh = cell(2,2);
    aod_directed  = cell(1,2);
    for i = 1:length(caseList)
        close;
        caseName = caseList{i}(1:end-4);
      %  targetpath = [root,'Input_0808_all\',caseName,'.nii'];

        vfxpath = [root,'Output\',suffix{i},caseName,'-_VelocityField_X.nii'];
        vfypath = [root,'Output\',suffix{i},caseName,'-_VelocityField_Y.nii'];
        vfzpath = [root,'Output\',suffix{i},caseName,'-_VelocityField_Z.nii'];
        aodpath = [root,'Output\',suffix{i},caseName,'-_TotalAOD.nii'];

        source = read_nii(sourcepath{i});
      %  target = read_nii(targetpath);
        % velocity field x
        vfx_ = read_nii(vfxpath);
        vfy_ = read_nii(vfypath);
        vfz_ = read_nii(vfzpath);
        scale_aod = read_nii(aodpath);
        
        vfx = vfx_(:,:,:,end);
        vfy = vfy_(:,:,:,end);
        vfz = vfz_(:,:,:,end);
        
        scale_aod = sqrt(vfx .^2 + vfy .^2 + vfz .^2);
        scale_aod(source==0)=0;
        load cover_height 
        load cover_area
 %       scale_aod = delete_cover(scale_aod,cover_height,cover_area);
        
        vectorx = vfx;
        vectory = vfy;
        vectorz = vfz;

        mesh{1,i} = perform_deform(source,vectorx,vectory,vectorz);
        
        direction(:,:,:,1) = vectorx;
        direction(:,:,:,2) = vectory;
        direction(:,:,:,3) = vectorz;

        pointnum = 50;
        centerdirection = zeros(size(direction));
        A = (1:pointnum) -1 ;
        B = (1:pointnum) -1 ;
        C = (1:pointnum) -1 ;

        [a,b,c] = ndgrid(A, B, C);
        pos = [a(:) b(:) c(:)];
        inorout = - pos + repmat([25 25 25], size(pos,1), 1);

        for a = 1:pointnum
            for b = 1:pointnum
              centerdirection(:,b,a,:) = inorout(1+pointnum*(b-1)+pointnum*pointnum*(a-1):pointnum*b+pointnum^2*(a-1),:);
            end
        end
        
        dir_ = centerdirection .* direction;
        direction = dir_(:,:,:,1) + dir_(:,:,:,2) + dir_(:,:,:,3);
        ddd = zeros(size(direction));
        ddd (direction  < 0) = 1; 
        ddd (direction  > 0) = -1;
        cur_aod = scale_aod .* ddd;
 
 %       direction(direction < 0) = 1;
 %       direction(direction > 0) = -1;
 %       cur_aod = scale_aod .* direction;
 
        cur_aod = to_plot_point_cloud(cur_aod, -1e5);
        idx1 = find(cur_aod(:,4)>1e-3);idx2 = find(cur_aod(:,4)<-1e-3);
        aod_directed{i} = cur_aod(union(idx1,idx2), :);
    end

    endo = mesh{1,1}; epi = mesh{1,2}; 
    aod1 = aod_directed{1}; aod2 = aod_directed{2};
    [endo, epi, aod1, aod2] = adjust_position(endo, epi, aod1, aod2);
    [class_endo, class_epi, class_aod1, class_aod2] = get_17_area_from_mesh(endo, epi, aod1, aod2);
    
 %   save([respath,'all_vari']);
    % close;
    figure(1458)
    hold on;
    scatter3(aod1(:,1),aod1(:,2),aod1(:,3));
    scatter3(aod2(:,1),aod2(:,2),aod2(:,3));
   % scatter3(endo(:,1),endo(:,2),endo(:,3));
   % scatter3(epi(:,1),epi(:,2),epi(:,3));
    %patch(endo,'FaceColor',[0.6 0.6 0.6],'facecolor','cyan','FaceAlpha',0.2,'EdgeAlpha',0.3);
    %patch(epi,'FaceColor',[0.6 0.6 0.6],'facecolor','cyan','FaceAlpha',0.2,'EdgeAlpha',0.3);
  %  pause()
end
