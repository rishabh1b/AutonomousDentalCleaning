x = 200;
y = 100;
rad = 10;
[xgrid, ygrid] = meshgrid(1:size(im,2), 1:size(im,1));
mask = ((xgrid-x).^2 + (ygrid-y).^2) <= rad.^2;
figure
imshow(mask)
values = im(mask);