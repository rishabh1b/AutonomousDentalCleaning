%% Read the Imge
im = imread('image2.JPG');
im_gray = rgb2gray(im);

%imtool(im_gray);
%% Put the corner points of transparent acrylic sheet as seen in the image
% square_clock = [202, 314;...
%                 606, 297;...
%                 657, 698;...
%                 224, 736];

% square_clock = [1102, 888;...
%                 2498, 607;...
%                 2695, 2197;...
%                 1193, 2310];
            
square_clock = [387, 320;...
                998, 293;...
                1134, 902;...
                341, 940];
            
warp_image_sz = 510;
front_facing_marker = [10 10; warp_image_sz 10; warp_image_sz warp_image_sz; 10 warp_image_sz];
homo = homography(square_clock, front_facing_marker);
image_size = size(im_gray);

%% Find the correspondences for each interior point
interior_pts = calculate_interior_pts(image_size,square_clock);

warped_points = warp_points(homo, interior_pts);

warped_points = ceil(warped_points);

%% Replace the pixel value in the rectified image
curr_image = im2double(im_gray);
rect_image = zeros(warp_image_sz+10,warp_image_sz+10);
ind_rect_image = sub2ind([warp_image_sz+10,warp_image_sz+10],warped_points(:,2),warped_points(:,1));
ind_vid_image = sub2ind(image_size, interior_pts(:,2), interior_pts(:,1));

rect_image(ind_rect_image) = curr_image(ind_vid_image);

imshow(rect_image);
%% Get the cropped version
[imcropped, rect] = imcrop(rect_image);
%%
figure
imshow(imcropped)
% imcrop_2 = imrotate(imcropped,90);
% figure
% imshow(imcrop_2)
%%
% Call RealTeethGetGumlinePoints2.m after this
RealTeethGetGumlinePoints2