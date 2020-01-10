function output = perform_deform(source,vectorx,vectory,vectorz)
    
    threshold_deform = -1e5;
    pc = to_plot_point_cloud(source,threshold_deform);
    %pc = source;
    movex = to_plot_point_cloud(vectorx, threshold_deform);
    movey = to_plot_point_cloud(vectory, threshold_deform);
    movez = to_plot_point_cloud(vectorz, threshold_deform);
    deformed = zeros(size(pc));
    
    deformed(:,1:3) = pc(:,1:3) + [movex(:,4) movey(:,4) movez(:,4)];
    deformed(:,4) = pc(:,4);
    
    pc = pc(~(pc(:,4)==0),:);
    deformedpc = deformed(~(deformed(:,4)==0),:);
    
    
    mesh.faces = boundary(deformedpc(:,1),deformedpc(:,2),deformedpc(:,3),0.65);
    tmp = repmat(deformedpc(1,1:3),length(deformedpc),1);
    tmp(unique(mesh.faces),:) = deformedpc(unique(mesh.faces),1:3);
    mesh.vertices = tmp;
 %   mesh.vertices = lpflow_trismooth(mesh.vertices,mesh.faces); 
    output = tmp;
end