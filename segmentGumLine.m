    for i = 46:1:90
        filename = strcat('dried_teeth_frontal\proc_',sprintf('%d',i),'.tiff');
        im = imread(filename);
        B = imsharpen(im);
        C = imadjust(B, stretchlim(B));
        %imtool(C)
        im_hsv = rgb2hsv(C);
        im_h = im_hsv(:,:,1);
        im_s = im_hsv(:,:,2);
        im_v = im_hsv(:,:,3);
        [m,n] = size(im_h);
        %% Thresholding in HSV space
        im_b = im_h > 0.8 & im_s > 0.2;
        im_b_2 = imfill(im_b, 'holes');
        %% Morphological cleaning
        struct = strel('disk',3);
        im_b_3 = imclose(im_b_2,struct);
        im_b_4 = imfill(im_b_3, 'holes');
        im_b_5 = imopen(im_b_4, struct);
        im_b_6 = edge(im_b_5);
        dilatedImage = imdilate(im_b_6,strel('disk',7));
        thinedImage = bwmorph(dilatedImage,'thin',inf);
        %imtool(thinedImage)
        %% Save the Output
    %     newFile = fullfile('outputs',filename);
    %     imwrite(thinedImage, newFile, 'tiff');
        %imtool(im_hsv(:,:,1))
        %imtool(im_hsv(:,:,2))
        %imtool(im_hsv(:,:,3))

        %% Get Only the gum Profile
        CC = bwconncomp(thinedImage);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [~, ind] = sort(numPixels, 'descend');
        %[biggest,idx] = max(numPixelsSorted);
        newImage = thinedImage;
        newImage(:) = 0;
        newImage(CC.PixelIdxList{ind(1)}) = 1;
        newImage(CC.PixelIdxList{ind(2)}) = 1;
        % TODO - If any pixel from top two are on both sides of the
        % centreline take the third one discarding the one closer to the
        % centreline
        %% Overlay on Original Image?
        CC = bwconncomp(dilatedImage);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [numPixelsSorted, ind] = sort(numPixels, 'descend');
        %[biggest,idx] = max(numPixelsSorted);
        newImageDilate = dilatedImage;
        newImageDilate(:) = 0;
        newImageDilate(CC.PixelIdxList{ind(1)}) = 1;
        newImageDilate(CC.PixelIdxList{ind(2)}) = 1;
        imgray = rgb2gray(im);
        Ifilt = imfilter(imgray,fspecial('gaussian'));
        im(newImageDilate) = Ifilt(newImageDilate);
        %% Show results side by side
        figure(1)
        subplot(1,2,1)
        imshow(im)
        subplot(1,2,2)
        imshow(newImage)
        %hgexport(gcf, fullfile('outputs2',filename), hgexport('factorystyle'), 'Format', 'jpeg');
    end