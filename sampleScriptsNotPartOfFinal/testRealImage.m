img = imread('test.jpg');
im_cropped = imcrop(img);
%%
imwrite(im_cropped, 'testcropped.jpg');

