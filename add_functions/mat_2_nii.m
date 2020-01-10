function pointcloud = mat_2_nii(mesh)

    pointnum = 50;unit=1;
    xmin = 0;ymin=0;zmin=0;
    pc = zeros(pointnum,pointnum,pointnum);
    A = xmin+((1:pointnum) -1 )*unit;
    B = ymin+((1:pointnum) -1 )*unit;
    C = zmin+((1:pointnum) -1 )*unit;

    [a,b,c] = ndgrid(A, B, C);
    pos = [a(:) b(:) c(:)];
    inorout = inpolyhedron(mesh,pos);
    
    for a = 1:pointnum
        for b = 1:pointnum
          pointcloud(:,b,a) = inorout(1+pointnum*(b-1)+pointnum*pointnum*(a-1):pointnum*b+pointnum^2*(a-1))*255;
        end
    end
end