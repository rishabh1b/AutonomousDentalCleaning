%% Get the tooth profile
show_figure_step = true;
threshold_num_pixels = 2000; %Avoid blobs smaller than this value
                             % Subject to change based on our experimental
                             % setup
margin_pixels_cropping = 0;
resolution = 5;
base_path = 'mseroutputs\img';
global count;
for i = count:1:count %52
    filename = strcat('dried_teeth_frontal\proc_',sprintf('%d',i),'.tiff');
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
    if show_figure_step
        figure
        imshow(M_roi_holes_filled)
    end
    seD = strel('disk',1);
    BWfinal = imerode(M_roi_holes_filled,seD);
    BWfinal = imerode(BWfinal,seD);
    BWfinal = bwareaopen(BWfinal,1000);
    if show_figure_step
        figure
        imshow(BWfinal)
        hold on
    end
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
    % Hack
    BWfinal(min_row_window-20:max_row_window+20,:) = 0;
    if show_figure_step
        figure
        imshow(BWfinal)
        title('Window Applied')
    end
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
    figure(2), imshow(bw_2), hold on
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
if show_figure_step
    figure
    imshow(bw_2)
    title('Horizontal Pixels Removed')
end
%% Get the median of the y points
%y_med = round(median(y_coords));
y_med = m/2;
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
if show_figure_step
    figure
    imshow(bw_upper)
    figure
    imshow(bw_lower)
end
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
if show_figure_step
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

    figure
    imshow(C)
    hold on
    for k = 1 : sz_1
        plot(selected_cols_upper{k},selected_rows_upper{k}, 'ro');
    end

    for k = 1:sz_2
        plot(selected_cols_lower{k},selected_rows_lower{k}, 'ro');
    end
    hold off
end
%% Peter corke trial
% [label,n] = imser(D, 'light');
% idisp(label)
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
first_elems_cols = zeros(sz_2,1);
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


figure(3)
imshow(C)
hold on
% Plot the points
for k = 1 : sz_1
    temp_sz = size(selected_cols_upper_sorted_cleaned{k},2);
    for t = 1 : temp_sz
        plot(selected_cols_upper_sorted_cleaned{k}(t),selected_rows_upper_sorted_cleaned{k}(t), 'ro');
        %pause(0.2)
    end
end

for k = 1 : sz_2
    temp_sz = size(selected_cols_lower_sorted_cleaned{k},2);
    for t = 1 : temp_sz
        plot(selected_cols_lower_sorted_cleaned{k}(t),selected_rows_lower_sorted_cleaned{k}(t), 'ro');
        %pause(0.2)
    end
end


% for k = 1 : sz_1
%     plot(selected_cols_upper_sorted{k},selected_rows_upper_sorted{k}, 'ro');
%% Get All the points for laser coverage
gumLinePointsUpper = [];
gumLinepointsLower = [];
for k = 1 : sz_1
   %gumLinePointsUpper = [gumLinePointsUpper;[selected_cols_upper_sorted_cleaned{k}(:),selected_rows_upper_sorted_cleaned{k}(:)]];
   gumLinePointsUpper = [gumLinePointsUpper;[selected_cols_upper{k}(:),selected_rows_upper{k}(:)]];
end
for k = 1 : sz_2
   %gumLinepointsLower = [gumLinepointsLower;[selected_cols_lower_sorted_cleaned{k}(:),selected_rows_lower_sorted_cleaned{k}(:)]];
   gumLinepointsLower = [gumLinepointsLower;[selected_cols_lower{k}(:),selected_rows_lower{k}(:)]];
end
end
