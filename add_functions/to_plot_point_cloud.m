function out = to_plot_point_cloud(res_plot, threshhold)

size_res = size(res_plot);

A = 1:size_res(1);
B = 1:size_res(2);
C = 1:size_res(3);

[a,b,c] = ndgrid(A, B, C);
sca = [a(:) b(:) c(:)];

ptsize = zeros(size(sca,1),1);
ptsizetop = zeros(size(sca,1),1);

for a = 1:size(res_plot,1)
    for b = 1:size(res_plot,2)
      ptsize(1+50*(b-1)+50*50*(a-1):50*b+50*50*(a-1)) = res_plot(:,b,a);
    end
end

sca_nozero = sca(ptsize > threshhold,:);
size_nozero = ptsize(ptsize > threshhold);

out = [sca_nozero, size_nozero];


end