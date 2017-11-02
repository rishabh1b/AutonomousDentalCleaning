show_figure_step = true;
threshold_num_pixels = 2000; %Avoid blobs smaller than this value
                             % Subject to change based on our experimental
                             % setup
margin_pixels_cropping = 0;
resolution = 30;
base_path = 'mseroutputs\img';
for i = 47:1:47 %52
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
    figure
    imshow(BWfinal)
    %%
      im_s = im_hsv(:,:,2);
      im_h = im_hsv(:,:,1);
%     im_b = im(:,:,3);
%     %im_s_c = imadjust(im_s);
     imtool(im_s)
     imtool(im_h)
%     imtool(im_h)
    imtool(C)
end
