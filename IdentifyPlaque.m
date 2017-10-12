show_figure_step = false;
threshold_num_pixels = 2000; %Avoid blobs smaller than this value
                             % Subject to change based on our experimental
                             % setup
margin_pixels_cropping = 1;
for i = 48:1:48 %52
    filename = strcat('dried_teeth_frontal\proc_',sprintf('%d',i),'.tiff');
    im = imread(filename);
    B = imsharpen(im);
    C = imadjust(B, stretchlim(B));
    %imtool(C)
    if show_figure_step
        figure,
        imshow(C);
    end
    %% Thresholding in HSV space
    %     im_hsv = rgb2hsv(C);
    %     im_h = im_hsv(:,:,1);
    %     im_s = im_hsv(:,:,2);
    %     im_v = im_hsv(:,:,3);
    %     [m,n] = size(im_h);
    %     im_b = im_h > 0.8 & im_s > 0.2;
    %     im_b_2 = imfill(im_b, 'holes');
    %% Morphological cleaning
    %     struct = strel('disk',8);
    %     im_b_3 = imclose(im_b_2,struct);
    %     im_b_4 = imfill(im_b_3, 'holes');
    %     im_b_5 = imopen(im_b_4, struct);
    %     %         imshow(im_b_5);
    %     for i =1:m
    %         for j=1:n
    %             if (im_b_5(i,j)==1)
    %                 C(i,j)=1;
    %             end
    %         end
    %     end
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
    %% Get rid of unwanted blobs
    BWfinal = bwareaopen(BWfinal, threshold_num_pixels);
%     label_w = bwlabel(BWfinal);
%     stats_w = regionprops(logical(BWfinal), 'Area');
%     areas = zeros(length(stats_w),1);
%     for j = 1 : length(stats_w)
%         areas(j) = stats_w(j).Area;
%     end
%     [~,ind] = min(areas);
%     BWfinal(label_w == ind) = 0;
    
    if show_figure_step
        figure
        imshow(BWfinal)
    end
    [row, ~] = find(BWfinal);
    row_min = min(row);
    row_max = max(row);
    
    row_min = max(row_min-margin_pixels_cropping, 1);
    row_max = min(row_max+margin_pixels_cropping,m);
    
    
    im_only_interest_area = C;
    im_only_interest_area(1:row_min,:) = 0;
    im_only_interest_area(row_max:m,:) = 0;
    im_r = im_only_interest_area(:,:,1);
    im_g = im_only_interest_area(:,:,2);
    im_b = im_only_interest_area(:,:,3);
    im_r(BWfinal == 1) = 0;
    im_g(BWfinal == 1) = 0;
    im_b(BWfinal == 1) = 0;
    im_segmented_plaque = cat(3,im_r,im_g, im_b);
       if show_figure_step
           figure
           imshow(im_segmented_plaque)
       end
    im_segmented_plaque_gray = rgb2gray(im_segmented_plaque);   
%% TODO Use Histogram function instead by grouping grayscale values in a vector
   segmented_plaque_1d = reshape(im_segmented_plaque_gray, [numel(im_segmented_plaque_gray),1]);
   [~,~,v] = find(segmented_plaque_1d);
   figure
   histogram(v);
%    thresh = graythresh(im_segmented_plaque_gray);
%    im_bw = imbinarize(im_segmented_plaque_gray,thresh);
%    figure
%    imshow(im_bw)
%% TODO get the hue and sat plane and look at the histogram plots for these
im_h_pl = im_segmented_plaque(:,:,1);
im_s_pl = im_segmented_plaque(:,:,2);

figure
imshow(im_segmented_plaque)
figure
imshow(im_h_pl)
figure
imshow(im_s_pl)
end