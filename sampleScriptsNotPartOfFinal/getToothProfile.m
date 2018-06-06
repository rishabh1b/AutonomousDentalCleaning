show_figure_step = false;
threshold_num_pixels = 2000; %Avoid blobs smaller than this value
                             % Subject to change based on our experimental
                             % setup
margin_pixels_cropping = 0;
resolution = 30;
for i = 49:1:49 %52
    filename = strcat('dried_teeth_frontal\proc_',sprintf('%d',i),'.tiff');
    im = imread(filename);
    B = imsharpen(im);
    C = imadjust(B, stretchlim(B));
    %imtool(C)
    if show_figure_step
        figure,
        imshow(C);
    end
    %% Find MSER regions
    D=rgb2gray(C);
    [m,n] = size(D);
    delta = 8;
    maxarea = 0.5;
    minarea = 0.0001;
    [r,~] = vl_mser(D,'MinDiversity',0.7,...
        'MaxVariation',0.2,...
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
    %M(1:m,:) = M(1:m,:) + M_roi;
    if show_figure_step
        %figure,imshow(M_roi);
    end
    %% Morphological Cleaning on the ROI image
    M_roi_holes_filled = imfill(M_roi,'holes');
    if show_figure_step
        figure, imshow(M_roi_holes_filled)
    end
    seD = strel('disk',1);
    BWfinal = imerode(M_roi_holes_filled,seD);
    BWfinal = imerode(BWfinal,seD);
    if show_figure_step
        figure, imshow(BWfinal), title('segmented image');
    end
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
end
%% Get the median of the y points
y_med = round(median(y_coords));
bw_2(y_med-40:y_med+40,:) = 0;
imshow(bw_2)
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
imshow(C)
hold on
for k = 1 : sz_1
    plot(selected_cols_upper{k},selected_rows_upper{k}, 'ro');
end

for k = 1:sz_2
    plot(selected_cols_lower{k},selected_rows_lower{k}, 'bo');
end
hold off
end

%%