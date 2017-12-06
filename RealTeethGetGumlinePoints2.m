% This script depends on RealTeethGetGumlinePoints.m.
% The script will take the cropped image (as obtained after homography) and
% give upper teeth points(only upper) in ascending order for now

resolution = 30;
B = imsharpen(imcropped);
C = B;
%imtool(C)
%% Find MSER regionss
D=im2uint8(C);
[m,n] = size(D);
delta = 3;
maxarea = 0.8;
minarea = 0.00001;
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
labels = bwlabel(M_roi);
stats = regionprops(M_roi,'Centroid');
for k = 1 : length(stats)
    if stats(k).Centroid(2) < 0.2 * m || stats(k).Centroid(2) > 0.8 * m
        M_roi(labels == k) = 0;
    end
end
figure,imshow(M_roi)
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
BWfinal = bwareaopen(BWfinal,300);
% figure
% imshow(BWfinal)
% hold on
%% Get the Rectangular window to separate out the upper and lower teeth
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
bw_2 = edge(BWfinal, 'Canny');
imshow(bw_2)
%% Hough Transform
[H,T,R] = hough(bw_2,'Theta',-90:0.5:-70);
P  = houghpeaks(H,45,'threshold',ceil(0.3*max(H(:))));
lines = houghlines(bw_2,T,R,P,'FillGap',2,'MinLength',2);
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
% figure
% imshow(bw_2)
% title('Horizontal Pixels Removed')
%% Get the median of the y points
y_med = round(median(y_coords));
%bw_2(min_row_window:max_row_window,:) = 0;
% bw_2(y_med-40:y_med+40,:) = 0;
% imshow(bw_2)
%% Get the points on the tooth profile
label_w = bwlabel(bw_2);
stats = regionprops(bw_2,'Centroid');
% First Sort into upper teeth and lower teeth
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
% figure
% imshow(bw_upper)
% figure
% imshow(bw_lower)
%% Get hold of points at a resolution in x-deirection along the teeth profile and plot
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
    [coS,ind] = sort(col);
    rwS = row(ind);
    numPoints = numel(coS);
    coS_T = coS';
    rwS_T = rwS';
    selected_cols_lower{k} = coS_T(1 : resolution: numPoints);
    selected_rows_lower{k} = rwS_T(1 : resolution: numPoints);
end
%% Plot
figure
imshow(bw_2)
title('Upper and Lower Teeth with Different Colours')
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
title('Upper and Lower Teeth with Different Colours')
hold on
for k = 1 : sz_1
    plot(selected_cols_upper{k},selected_rows_upper{k}, 'ro');
end

for k = 1:sz_2
    plot(selected_cols_lower{k},selected_rows_lower{k}, 'bo');
end
hold off

% Trying to get the points in order
%% Get the Points on the upper teeth in order
selected_cols_upper_sorted = cell(sz_1,1);
selected_rows_upper_sorted = cell(sz_1,1);
selected_cols_lower_sorted = cell(sz_2,1);
selected_rows_lower_sorted = cell(sz_2,1);

% Sort the First element in each column of upper teeth
first_elems_cols = zeros(sz_1,1);
for k = 1 : sz_1
    first_elems_cols(k) = selected_cols_upper{k}(1);
end
[~, ind] = sort(first_elems_cols);
for k = 1 : sz_1
    curr_cols = selected_cols_upper{ind(k)};
    curr_rows = selected_rows_upper{ind(k)};
    [sorted_cols,ind2] = sort(curr_cols);
    selected_cols_upper_sorted{k} = sorted_cols;
    selected_rows_upper_sorted{k} = curr_rows(ind2);
end

% Sort the First column in each column of lower teeth
for k = 1 : sz_2
    first_elems_cols(k) = selected_cols_lower{k}(1);
end
[~, ind] = sort(first_elems_cols);
for k = 1 : sz_2
    curr_cols = selected_cols_lower{ind(k)};
    curr_rows = selected_rows_lower{ind(k)};
    [sorted_cols,ind2] = sort(curr_cols);
    selected_cols_lower_sorted{k} = sorted_cols;
    selected_rows_lower_sorted{k} = curr_rows(ind2);
end

%% Clean the Sorted Rows and Columns to get nice ascending points
selected_cols_upper_sorted_cleaned = cell(sz_1,1);
selected_rows_upper_sorted_cleaned = cell(sz_1,1);
selected_cols_lower_sorted_cleaned = cell(sz_2,1);
selected_rows_lower_sorted_cleaned = cell(sz_2,1);

min_last_tooth = 0;
max_last_tooth = 0;
% First clean the upper teeth
for k = 1:sz_1
    rows = selected_rows_upper_sorted{k};
    cols = selected_cols_upper_sorted{k};
    if cols(1) > min_last_tooth && cols(1) < max_last_tooth
        continue;
    end
    curr_sz = size(rows,2);
    cleaned_rows = rows(1);
    cleaned_cols = cols(1);
    last_added_row_elem = rows(1);
    [~, ind] = min(rows);
    for j = 1 : ind-1
        if rows(j+1) > last_added_row_elem || abs(rows(j+1) - last_added_row_elem) > 80
            continue;
        else
            cleaned_rows = [cleaned_rows rows(j+1)];
            cleaned_cols = [cleaned_cols cols(j+1)];
            last_added_row_elem = rows(j+1);
        end
    end
    cleaned_rows = [cleaned_rows rows(ind)];
    cleaned_cols = [cleaned_cols cols(ind)];
    last_added_row_elem = rows(ind);
    for j = ind : curr_sz - 1
        if rows(j+1) < last_added_row_elem || abs(rows(j+1) - last_added_row_elem) > 80
            continue;
        else
            cleaned_rows = [cleaned_rows rows(j+1)];
            cleaned_cols = [cleaned_cols cols(j+1)];
            last_added_row_elem = rows(j+1);
        end
    end
selected_cols_upper_sorted_cleaned{k} = cleaned_cols;
selected_rows_upper_sorted_cleaned{k} = cleaned_rows;
min_last_tooth = min(cleaned_cols);
max_last_tooth = max(cleaned_cols);
end

min_last_tooth = 0;
max_last_tooth = 0;

% Then clean the lower teeth
for k = 1:sz_2
    rows = selected_rows_lower_sorted{k};
    cols = selected_cols_lower_sorted{k};
    if cols(1) > min_last_tooth && cols(1) < max_last_tooth
        continue;
    end
    curr_sz = size(rows,2);
    cleaned_rows = rows(1);
    cleaned_cols = cols(1);
    last_added_row_elem = rows(1);
    [~, ind] = sort(rows, 'descend');
    if numel(ind) > 1
        ind = ind(2);
    end

    for j = 1 : ind-1
        if rows(j+1) < last_added_row_elem || abs(rows(j+1) - last_added_row_elem) > 200
            continue;
        else
            cleaned_rows = [cleaned_rows rows(j+1)];
            cleaned_cols = [cleaned_cols cols(j+1)];
            last_added_row_elem = rows(j+1);
        end
    end
    cleaned_rows = [cleaned_rows rows(ind)];
    cleaned_cols = [cleaned_cols cols(ind)];
    last_added_row_elem = rows(ind);
    for j = ind : curr_sz - 1
        if rows(j+1) > last_added_row_elem || abs(rows(j+1) - last_added_row_elem) > 100
            continue;
        else
            cleaned_rows = [cleaned_rows rows(j+1)];
            cleaned_cols = [cleaned_cols cols(j+1)];
            last_added_row_elem = rows(j+1);
        end
    end
selected_cols_lower_sorted_cleaned{k} = cleaned_cols;
selected_rows_lower_sorted_cleaned{k} = cleaned_rows;
min_last_tooth = min(cleaned_cols);
max_last_tooth = max(cleaned_cols);
end

% Animate the Points in ascending order in the cropped image
% figure
% imshow(C)
% hold on
% % Plot the points
% for k = 1 : sz_1
%     temp_sz = size(selected_cols_upper_sorted_cleaned{k},2);
%     for t = 1 : temp_sz
%         plot(selected_cols_upper_sorted_cleaned{k}(t),selected_rows_upper_sorted_cleaned{k}(t), 'ro');
%         pause(0.2)
%     end
% end
% 
% for k = 1 : sz_2
%     temp_sz = size(selected_cols_lower_sorted_cleaned{k},2);
%     for t = 1 : temp_sz
%         plot(selected_cols_lower_sorted_cleaned{k}(t),selected_rows_lower_sorted_cleaned{k}(t), 'go');
%         pause(0.2)
%     end
% end

%% Get just the Upper Teeth Points and scale it back to original image
upper_teeth_points = [];
for k = 1 : sz_1
    temp_sz = size(selected_cols_upper_sorted_cleaned{k},2);
    for t = 1 : temp_sz
        upper_teeth_points = [upper_teeth_points;[selected_cols_upper_sorted_cleaned{k}(t),selected_rows_upper_sorted_cleaned{k}(t)]];
    end
end

upper_teeth_points(:,1) = upper_teeth_points(:,1) + rect(1);
upper_teeth_points(:,2) = upper_teeth_points(:,2) + rect(2);
%% Plot and confirm in Front Facing Image
figure
imshow(rect_image)
title('Points on the gumline of upper tooth in rectified image')
hold on
plot(upper_teeth_points(:,1), upper_teeth_points(:,2), 'r*')

%% Do Inverse Homography to transfer the points on the original image
[m,~] = size(upper_teeth_points);
upper_teeth_points_h = [upper_teeth_points';ones(1,m)];
upper_teeth_points_h_orig = homo \ upper_teeth_points_h;
upper_teeth_points_h_orig = upper_teeth_points_h_orig ./(upper_teeth_points_h_orig(3,:));
upper_teeth_points_orig = upper_teeth_points_h_orig(1:2,:);
upper_teeth_points_orig = upper_teeth_points_orig';
figure
imshow(im)
title('Points on the gumline of upper tooth in Original image')
hold on
for i = 1 : m
    plot(upper_teeth_points_orig(i,1),upper_teeth_points_orig(i,2),'r*')
    pause(0.2)
end