show_figure_step = true;
for i = 48:1:48
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
        figure,imshow(M_roi);
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
end