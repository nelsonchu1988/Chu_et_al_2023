function istogramma(x, y, b, x_max, x_min, scale)

%x = bm.r(:,1);
%y = bm.r(:,2);

dx = (x_max-x_min)/b;

xq = ceil((x-x_min)/dx);
yq = ceil((y-x_min)/dx);

p_xy = zeros(b,b);

for i = min(xq):1:max(xq)
    xi = find(xq==i);
    for j = min(yq):1:max(yq)
        p_xy(i,j) = numel(find(yq(xi)==j));
    end
end

xp = scale*(x_min:dx:x_max-dx);
surf(xp, xp, -p_xy)
colormap gray
view(2)
shading flat
xlim(scale*[x_min x_max])
ylim(scale*[x_min x_max])

axis off