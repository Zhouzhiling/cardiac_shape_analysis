function mesh = deformDiscrepancy(mesh, deformSize)
    pointNum = length(mesh.vertices(:,1));
    %R = randi(7);
    R = 3;
    if(~(R == 2))
        for i = 1:pointNum;
            % towards the point
            mesh.vertices(i,:) = mesh.vertices(i,:) + deformSize .* mesh.diffVector(i,:); 
            % in the normal direction
            %mesh.vertices(i,:) = mesh.vertices(i,:) + deformSize .* mesh.deformScale(i) .* mesh.normal(i,:); 
        end
    else
        for i = 1:pointNum;
            % towards the point
            %mesh.vertices(i,:) = mesh.vertices(i,:) + deformSize .* mesh.diffVector(i,:); 
            % in the normal direction
            mesh.vertices(i,:) = mesh.vertices(i,:) + deformSize .* mesh.deformScale(i) .* mesh.normal(i,:); 
        end
    end
end