%% from vtk to nii
% 从vtk生成nii文件
root = 'C:\Users\37908\Desktop\tmp\out\AWS_OtherCases_Average\';
outroot = 'C:\Users\37908\Desktop\tmp\input_data_nii_shell\EpiED\';
filename = 'DeterministicAtlas__EstimatedParameters__Template_othercase.vtk';

[vertex, face] = read_vtk(strcat(root,filename));
mesh.vertices = vertex';
mesh.faces = face';

pointsize = 50;
xyzSize = [92.7332   94.7080  108.8858];
minPos = [min(mesh.vertices(:,1)),min(mesh.vertices(:,2)),min(mesh.vertices(:,3))];
mesh.vertices = mesh.vertices - repmat(minPos, length(mesh.vertices(:,1)),1);

maxLength = max(xyzSize);
unitDist = maxLength / pointsize;
pointcloud = zeros(pointsize,pointsize,pointsize);
for i = 1:pointsize
    i
    for j = 1:pointsize
        j
        for k = 1:pointsize
            pos = [i*unitDist, j*unitDist, k*unitDist];
            if(inpolyhedron(mesh,pos))
                pointcloud(i,j,k) = 255;
            end
        end
    end
end
y=make_nii(pointcloud);
save_nii(y,strcat(outroot,filename,'.nii'));

