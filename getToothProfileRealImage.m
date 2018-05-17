%% Get the tooth profile
show_figure_step = false;
threshold_num_pixels = 2000; %Avoid blobs smaller than this value
                             % Subject to change based on our experimental
                             % setup
margin_pixels_cropping = 0;
resolution = 30;
base_path = 'mseroutputs\img';
for i = 49:1:49 %52
    %filename = strcat('dried_teeth_frontal\proc_',sprintf('%d',i),'.tiff');
    filename = 'testcropped.jpg';
    im = imread(filename);
    B = imsharpen(im);
    %C = imadjust(B, stretchlim(B));
    C = B;
    im_hsv = rgb2hsv(C);
    %imtool(C)
    if show_figure_step
        figure,
        imshow(C);
    end
    %% Find MSER regionss
    D=rgb2gray(C);
    [m,n] = size(D);
    delta = 8;
    maxarea = 0.5;
    minarea = 0.0001;
    [r,~] = vl_mser(D,'MinDiversity',0.7,...
        'MaxVariation',0.4,...
        'Delta',delta, 'DarkOnBright', 0, 'MaxArea', maxarea, 'MinArea', minarea ) ;
    sAll = [];
    for x=r'
        s = vl_erfill(D,x) ;
        sAll = [sAll;s];
    end
    % Obtain the output in the original size image
    % Use this if want to constraint MSER to only a subset of entire image
    %M = zeros(m,n) ;
    %     M_roi = false(m,n);
    %     M_roi(sAll) = 1;
    M_roi = false(m,n);
    M_roi(sAll) = 1;
    M_roi = bwareaopen(M_roi,1000);
    %M(1:m,:) = M(1:m,:) + M_roi;
    if show_figure_step
        %figure,imshow(M_roi);
    end
    labels = bwlabel(M_roi);
    stats = regionprops(M_roi,'Centroid');
    for k = 1 : length(stats)
        if stats(k).Centroid(2) < 0.2 * m || stats(k).Centroid(2) > 0.8 * m
            M_roi(labels == k) = 0;
        end
    end
%     imgname = strcat(base_path, sprintf('%d',i));
%     imgname = strcat(imgname, '.jpg');
%     imwrite(M_roi,imgname);
%     figure
%     imshow(M_roi)
    M_roi_holes_filled = imfill(M_roi,'holes');
    figure
    imshow(M_roi_holes_filled)
    seD = strel('disk',1);
    BWfinal = imerode(M_roi_holes_filled,seD);
    BWfinal = imerode(BWfinal,seD);
    BWfinal = bwareaopen(BWfinal,1000);
    figure
    imshow(BWfinal)
    hold on
    %% Get the Rectangular window
    labels = bwlabel(BWfinal);
    stats = regionprops(BWfinal, 'Centroid', 'BoundingBox');
    min_row_window = m;
    max_row_window = 1;
    for k = 1 : length(stats)
        [row,~] = find(labels == k);
        centre = stats(k).Centroid;
        bbox = stats(k).BoundingBox;
        if bbox(4) > m/2
            continue;
        end
        if stats(k).Centroid(2) <= m / 2
            plot(centre(1), centre(2),'r*')
            curr_max_row = max(row);
            if min_row_window > curr_max_row
                min_row_window = curr_max_row;
            end
        else
            plot(centre(1), centre(2),'g*')
            curr_min_row = min(row);
            if max_row_window < curr_min_row
                max_row_window = curr_min_row;
            end
        end
    end
    %% Apply the Rectangular window
    BWfinal(min_row_window:max_row_window,:) = 0;
    figure
    imshow(BWfinal)
    title('Window Applied')
    %%
      im_s = im_hsv(:,:,2);
      im_h = im_hsv(:,:,1);
%     im_b = im(:,:,3);
%     %im_s_c = imadjust(im_s);
     %imtool(im_s)
     %imtool(im_h)
%     imtool(im_h)
    %imtool(C)
    %%
     bw_2 = edge(BWfinal, 'Canny');
    %imshow(bw_2)
    %% Hough Transform
    [H,T,R] = hough(bw_2,'Theta',-90:0.5:-70);
    P  = houghpeaks(H,45,'threshold',ceil(0.3*max(H(:))));
    lines = houghlines(bw_2,T,R,P,'FillGap',20,'MinLength',25);
    figure, imshow(bw_2), hold on
max_len = 0;
y_coords = zeros(length(lines) * 2,1);
index = 1;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
   y_coords(index:index+1) = xy(:,2);
   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
   index = index + 2;
   bw_2(floor(xy(1,2))-10 : ceil(xy(2,2))+10, floor(xy(1,1)) : ceil(xy(2,1))) = 0;
end
figure
imshow(bw_2)
title('Horizontal Pixels Removed')
%% Get the median of the y points
y_med = round(median(y_coords));
%bw_2(min_row_window:max_row_window,:) = 0;
% bw_2(y_med-40:y_med+40,:) = 0;
% imshow(bw_2)
%% Get the points on the tooth profile
label_w = bwlabel(bw_2);
stats = regionprops(bw_2,'Centroid');
% Sort into upper teeth and lower teeth
centroids = zeros(length(stats),2);
upperTeethLabels = [];
lowerTeethLabels = [];
bw_upper = false(size(bw_2));
bw_lower = false(size(bw_2));
for k = 1:length(stats)
    if stats(k).Centroid(2) < y_med
        bw_upper(label_w == k) = 1;
        upperTeethLabels = [upperTeethLabels;k];
    else
        bw_lower(label_w == k) = 1;
        lowerTeethLabels = [lowerTeethLabels;k];
    end
end
figure
imshow(bw_upper)
figure
imshow(bw_lower)
%% Get hold of points at a resolution on teeth and plot
sz_1 = numel(upperTeethLabels);
selected_cols_upper = cell(sz_1,1);
selected_rows_upper = cell(sz_1,1);
for k = 1:sz_1
    [row, col] = find(label_w == upperTeethLabels(k));
    %rwS = sort(row);
    [coS,ind] = sort(col);
    rwS = row(ind);
    numPoints = numel(coS);
    %numPointsPicked = numPoints / resolution;
    coS_T = coS';
    rwS_T = rwS';
    selected_cols_upper{k} = coS_T(1 : resolution: numPoints);
    selected_rows_upper{k} = rwS_T(1 : resolution: numPoints);
end

%% Repitition for now refactor into a function later
sz_2 = numel(lowerTeethLabels);
selected_cols_lower = cell(sz_2,1);
selected_rows_lower = cell(sz_2,1);
for k = 1:sz_2
    [row, col] = find(label_w == lowerTeethLabels(k));
    %rwS = sort(row);
    [coS,ind] = sort(col);
    rwS = row(ind);
    numPoints = numel(coS);
    %numPointsPicked = numPoints / resolution;
    coS_T = coS';
    rwS_T = rwS';
    selected_cols_lower{k} = coS_T(1 : resolution: numPoints);
    selected_rows_lower{k} = rwS_T(1 : resolution: numPoints);
end
%% Plot
figure
imshow(bw_2)
hold on
for k = 1 : sz_1
    plot(selected_cols_upper{k},selected_rows_upper{k}, 'ro');
end

for k = 1:sz_2
    plot(selected_cols_lower{k},selected_rows_lower{k}, 'go');
end
hold off
%%
figure
imshow(C)
hold on
for k = 1 : sz_1
    plot(selected_cols_upper{k},selected_rows_upper{k}, 'ro');
end

for k = 1:sz_2
    plot(selected_cols_lower{k},selected_rows_lower{k}, 'bo');
end
hold off
%% Peter corke trial
% [label,n] = imser(D, 'light');
% idisp(label)
%% Get the Points on the upper teeth in order
selected_cols_upper_sorted = cell(sz_1,1);
selected_rows_upper_sorted = cell(sz_1,1);
selected_cols_lower_sorted = cell(sz_2,1);
selected_rows_lower_sorted = cell(sz_2,1);

% Sort the First element in each column of upper teeth
% first_elems_cols = zeros(sz_1,1);
% for k = 1 : sz_1
%     first_elems_cols(k) = selected_cols_upper{k}(1);
% end
% [~, ind] = sort(first_elems_cols);
% selected_cols_upper_sorted{1:sz_1} = selected_cols_upper{ind};
% selected_rows_upper_sorted{1:sz_1} = selected_rows_upper{ind};
% 
% % Sort the First column in each column of lower teeth
% for k = 1 : sz_2
%     first_elems_cols(k) = selected_cols_lower{k}(1);
% end
% [~, ind] = sort(first_elems_cols);
% selected_cols_lower_sorted{1:sz_2} = selected_cols_lower{ind};
% selected_rows_lower_sorted{1:sz_2} = selected_rows_upper{ind};

figure
imshow(C)
hold on
% Plot the points
for k = 1 : sz_1
    temp_sz = size(selected_cols_upper{k},2);
    for t = 1 : temp_sz
        plot(selected_cols_upper{k}(t),selected_rows_upper{k}(t), 'ro');
        pause(1)
    end
end

for k = 1 : sz_2
    plot(selected_cols_lower{k},selected_rows_lower{k}, 'go');
    pause(1)
end

% for k = 1 : sz_1
%     plot(selected_cols_upper_sorted{k},selected_rows_upper_sorted{k}, 'ro');
%     pause(0.001)
% end
% 
% for k = 1 : sz_2
%     plot(selected_cols_lower_sorted{k},selected_rows_lower_sorted{k}, 'ro');
%     pause(0.001)
% end
end
