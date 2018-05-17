% im = imread('trial_2.png');
% im=imcrop(im);
% t = @(x) x(:,1).*-pi/10;
% f = @(x) [x(:,1).*cos(t(x))+x(:,2).*sin(t(x)),-x(:,1).*sin(t(x))+x(:,2).*cos(t(x))];
% g = @(x, unused) f(x);
% tform = maketform('custom', 2, 2, [], g, []);
% IM = imtransform(im, tform, 'UData', [-1 1], 'VData', [-1 1], ...
%    'XData', [-1 1], 'YData', [-1 1]);
% figure(2),
% imshow(IM)
%%
%% Read the Imge
%im = imread('trial_2.png');
%im = imread('image2.JPG');
%im = imread('newSetupImages/setup_3.jpg');
%im = imread('image_raw_screenshot_09.02.2018.png');
%im_gray = rgb2gray(im);

im_gray = imread('newsetup_2.png');
im = im_gray;
figure
imshow(im_gray);
%% Resize to be consistent with the cv_bridge images
%im_gray = imresize(im_gray,[964,1296]);
%imtool(im_gray)
%% Put the corner points of transparent acrylic sheet as seen in the image
% square_clock = [202, 314;...
%                 606, 297;...
%                 657, 698;...
%                 224, 736];

% square_clock = [1102, 888;...
%                 2498, 607;...
%                 2695, 2197;...
%                 1193, 2310];
            
% square_clock = [387, 320;...
%                 998, 293;...
%                 1134, 902;...
%                 341, 940];
            
% square_clock = [330, 258;...
%                 1002, 232;...
%                 1152, 960;...
%                 246, 995];

% Nice Calibration
% square_clock = [620, 481;...
%                 996, 468;...
%                 1057, 843;...
%                 599, 843];

% Perfect Top View
% square_clock = [291, 71;...
%                 1098, 37;...
%                 1124, 879;...
%                 316, 888];
% Last working
% square_clock = [289, 390;...
%                 752, 340;...
%                 848, 854;...
%                 296, 932];
% From Windows
% square_clock = [218, 58;...
%                 1012, 81;...
%                 1109, 807;...
%                 161, 815];
% New Setup          
square_clock = [189, 335;...
                867, 268;...
                945, 841;...
                145, 907];
%warp_image_sz = 510;
warp_image_sz = 310;
front_facing_marker = [10 10; warp_image_sz 10; warp_image_sz warp_image_sz; 10 warp_image_sz];
homo = homography(square_clock, front_facing_marker);
image_size = size(im_gray);

%% Find the correspondences for each interior point
interior_pts = calculate_interior_pts(image_size,square_clock);

warped_points = warp_points(homo, interior_pts);

warped_points = ceil(warped_points);

%% Replace the pixel value in the rectified image
margin_size = 5;
curr_image = im2double(im_gray);
rect_image = zeros(warp_image_sz+margin_size,warp_image_sz+margin_size);
ind_rect_image = sub2ind([warp_image_sz+margin_size,warp_image_sz+margin_size],warped_points(:,2),warped_points(:,1));
ind_vid_image = sub2ind(image_size, interior_pts(:,2), interior_pts(:,1));

rect_image(ind_rect_image) = curr_image(ind_vid_image);
figure
imshow(rect_image);
%% Get the cropped version
[imcropped, rect] = imcrop(rect_image);
t = @(x) x(:,1).*-pi/10;
f = @(x) [x(:,1).*cos(t(x))+x(:,2).*sin(t(x)),-x(:,1).*sin(t(x))+x(:,2).*cos(t(x))];
g = @(x, unused) f(x);
tform = maketform('custom', 2, 2, [], g, []);
imcropped = imtransform(imcropped, tform, 'UData', [-1 1], 'VData', [-1 1], ...
   'XData', [-1 1], 'YData', [-1 1]);
% figure(2),
%imshow(IM)
%%
figure
imshow(imcropped)
% imcrop_2 = imrotate(imcropped,90);
% figure
% imshow(imcrop_2)
%%
% Call RealTeethGetGumlinePoints2.m after this
RealTeethGetGumlinePoints2
%RealTeethGetGumlinePoints3