%% Get the tooth profile
show_figure_step = false;
threshold_num_pixels = 2000; %Avoid blobs smaller than this value
                             % Subject to change based on our experimental
                             % setup
margin_pixels_cropping = 0;
resolution = 30;
base_path = 'mseroutputs\img';
for i = 52:1:52 %52
    filename = strcat('dried_teeth_frontal\proc_',sprintf('%d',i),'.tiff');
    im = imread(filename);
    B = imsharpen(im);
    C = imadjust(B, stretchlim(B));
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
end
