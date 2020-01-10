function [class_endo, class_epi, class_aod1, class_aod2] = get_17_area_from_mesh(varargin)
    
    Endo = varargin{1};
    Epi = varargin{2};
       
    if(length(Endo)==1)
        endo = Endo;
        epi = Epi;        
    else
        endo.vertices = Endo;
        epi.vertices = Epi;
    end
    
     if nargin == 2
        aod1 = endo.vertices;
        aod2 = epi.vertices;
    elseif nargin ==4
        aod1 = varargin{3};
        aod2 = varargin{4};
     end
    
    size_endo = length(endo.vertices);
    size_epi = length(epi.vertices);
    size_aod1 = length(aod1);
    size_aod2 = length(aod2);
        
    % move central to (0 0 0)
    %central = [0 0 0]; 
    centralepi = mean(epi.vertices);
    epi.vertices = epi.vertices - repmat(centralepi,length(epi.vertices),1);
    aod2(:,1:3) = aod2(:,1:3) - repmat(centralepi,length(aod2),1);
    
%    centralendo = mean(endo.vertices);
    endo.vertices = endo.vertices - repmat(centralepi,length(endo.vertices),1);
    aod1(:,1:3) = aod1(:,1:3) - repmat(centralepi, length(aod1),1);

    xmin = min(epi.vertices(:,1)); xmax = max(epi.vertices(:,1));
    ymin = min(epi.vertices(:,2)); ymax = max(epi.vertices(:,2));
    zmin = min(epi.vertices(:,3)); zmax = max(epi.vertices(:,3));
    xymin = min(xmin,ymin); xymax = max(xmax,ymax);

    % cut horizontally
    basal_height = max(epi.vertices(:,3));
    bottom_height = min(endo.vertices(:,3));
    rest_height = basal_height - bottom_height;
    each_part_height = rest_height / 3;
    top_height = basal_height - each_part_height;
    middle_height = basal_height - 2*each_part_height;
    % Analytic:
    % z = middle_height
    % z = top_height
    % z = bottom_height
        
    % clear allepi allendo;
    % divide the vertices in height
%    pc = [epi.vertices;endo.vertices];
    pc_withaod = [endo.vertices; epi.vertices; aod1(:,1:3); aod2(:,1:3)];
    
    class = zeros(size(pc_withaod,1),1);
    class(pc_withaod(:,3) < bottom_height) = 17;
    class(pc_withaod(:,3) > top_height) = -1;
    idx1 = find(pc_withaod(:,3) > middle_height);
    idx2 = find(pc_withaod(:,3) <= top_height);
    idx3 = find(pc_withaod(:,3) >= bottom_height);
    idx4 = find(pc_withaod(:,3) <= middle_height);
    class(intersect(idx1,idx2)) = -2;
    class(intersect(idx3,idx4)) = -3;
        
    % Apical divide
        class((pc_withaod(:,1) <= pc_withaod(:,2))&...
                     (pc_withaod(:,1) < -pc_withaod(:,2))&...
                     (class==-3))        = 14;
        class((pc_withaod(:,1) > -pc_withaod(:,2))&...
                     (pc_withaod(:,1) <= pc_withaod(:,2))&...
                     (class==-3))        = 15;          
        class((pc_withaod(:,1) >= pc_withaod(:,2))&...
                     (pc_withaod(:,1) > -pc_withaod(:,2))&...
                     (class==-3))        = 16;    
        class((pc_withaod(:,1) < -pc_withaod(:,2))&...
                     (pc_withaod(:,1) >= pc_withaod(:,2))&...
                     (class==-3))        = 13;
                 
    % Mid divide
        class((pc_withaod(:,2) > 0)&...
                     (pc_withaod(:,1) >= pc_withaod(:,2)/3^(1/2))&...
                     (class==-2))        = 10;
        class((pc_withaod(:,2) <= 0)&...
                     (pc_withaod(:,1) > -pc_withaod(:,2)/3^(1/2))&...
                     (class==-2))        = 11;
        class((pc_withaod(:,2) < 0)&...
                     (pc_withaod(:,1) <= pc_withaod(:,2)/3^(1/2))&...
                     (class==-2))        = 7;
        class((pc_withaod(:,2) >= 0)&...
                     (pc_withaod(:,1) < -pc_withaod(:,2)/3^(1/2))&...
                     (class==-2))        = 8;
        class((pc_withaod(:,1) >= -pc_withaod(:,2)/3^(1/2))&...
                     (pc_withaod(:,1) < pc_withaod(:,2)/3^(1/2))&...
                     (class==-2))        = 9;
        class((pc_withaod(:,1) > pc_withaod(:,2)/3^(1/2))&...
                     (pc_withaod(:,1) <= -pc_withaod(:,2)/3^(1/2))&...
                     (class==-2))        = 12;
                 
        % Basal divide
        class((pc_withaod(:,2) > 0)&...
                     (pc_withaod(:,1) >= pc_withaod(:,2)/3^(1/2))&...
                     (class==-1))        = 4;
        class((pc_withaod(:,2) <= 0)&...
                     (pc_withaod(:,1) > -pc_withaod(:,2)/3^(1/2))&...
                     (class==-1))        = 5;
        class((pc_withaod(:,2) < 0)&...
                     (pc_withaod(:,1) <= pc_withaod(:,2)/3^(1/2))&...
                     (class==-1))        = 1;
        class((pc_withaod(:,2) >= 0)&...
                     (pc_withaod(:,1) < -pc_withaod(:,2)/3^(1/2))&...
                     (class==-1))        = 2;
        class((pc_withaod(:,1) >= -pc_withaod(:,2)/3^(1/2))&...
                     (pc_withaod(:,1) < pc_withaod(:,2)/3^(1/2))&...
                     (class==-1))        = 3;
        class((pc_withaod(:,1) > pc_withaod(:,2)/3^(1/2))&...
                     (pc_withaod(:,1) <= -pc_withaod(:,2)/3^(1/2))&...
                     (class==-1))        = 6;
        
    pc_withaod = [endo.vertices; epi.vertices; aod1(:,1:3); aod2(:,1:3)];
    
    
    class_endo = class(1:size_endo);
    class_epi =  class(size_endo+1:size_endo+size_epi);
    class_aod1 = class(size_endo+size_epi+1:size_endo+size_epi+size_aod1);
    class_aod2 = class(end-size_aod2+1:end);
    
   

end