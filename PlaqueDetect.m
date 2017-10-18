for i = 48:1:85
    filename = strcat('proc_',sprintf('%d',i),'.tiff');
    im = imread(filename);
    B = imsharpen(im);
    C = imadjust(B, stretchlim(B));
    %imtool(C)
    im_hsv = rgb2hsv(C);
    im_h = im_hsv(:,:,1);
    im_s = im_hsv(:,:,2);
    im_v = im_hsv(:,:,3);
    [m,n] = size(im_h);
%     figure,
%     imshow(C);
    %% Thresholding in HSV space
    im_b = im_h > 0.8 & im_s > 0.2;
    im_b_2 = imfill(im_b, 'holes');
    %% Morphological cleaning
    struct = strel('disk',8);
    im_b_3 = imclose(im_b_2,struct);
    im_b_4 = imfill(im_b_3, 'holes');
    im_b_5 = imopen(im_b_4, struct);
    %         imshow(im_b_5);
    for i =1:m
        for j=1:n
            if (im_b_5(i,j)==1)
                C(i,j)=1;
            end
        end
    end
    D=rgb2gray(C);
    delta = 8;
    maxarea = 0.5;
    minarea = 0.0001;
    [r,~] = vl_mser(D,'MinDiversity',0.7,...
        'MaxVariation',0.2,...
        'Delta',delta, 'DarkOnBright', 0, 'MaxArea', maxarea, 'MinArea', minarea ) ;
    M = zeros(m,n) ;
    sAll = [];
    for x=r'
        s = vl_erfill(D,x) ;
        sAll = [sAll;s];
    end
    % Obtain the output in the original size image
    M_roi = zeros(m,n);
    M_roi(sAll) = 1;
    M(1:m,:) = M(1:m,:) + M_roi;
    figure,imshow(M);
    r=C(:,:,1);
    g=C(:,:,2);
    b=C(:,:,3);
    struct1=strel('line',1,30);
    im_b= r>215 & g>150 & b<180;
    im_b=imfill(im_b,'holes');
    %im_b=imopen(im_b,struct1);
%     figure
%     imshow(im_b);
    figure(1)
    subplot(1,2,1)
    imshow(im)
    subplot(1,2,2)
    imshow(im_b)
    hgexport(gcf, fullfile('plaqueoutput',filename), hgexport('factorystyle'), 'Format', 'jpeg');
end