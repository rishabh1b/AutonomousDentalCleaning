for i = 48:1:48
    filename = strcat('dried_teeth_frontal\proc_',sprintf('%d',i),'.tiff');
    tic
    im = imread(filename);
    B = imsharpen(im);
    C = imadjust(B, stretchlim(B));
    im_hsv = rgb2hsv(C);
    im_h = im_hsv(:,:,1);
    im_s = im_hsv(:,:,2);
    im_v = im_hsv(:,:,3);
    [m,n] = size(im_h);
    %% Thresholding in HSV space
    im_b = im_h > 0.8 & im_s > 0.2;
    im_b_2 = imfill(im_b, 'holes');
    %% Morphological cleaning
    struct = strel('disk',8);
    im_b_3 = imclose(im_b_2,struct);
    im_b_4 = imfill(im_b_3, 'holes');
    im_b_5 = imopen(im_b_4, struct);
    %         imshow(im_b_5);
    for k =1:m
        for j=1:n
            if (im_b_5(k,j)==1)
                C(k,j)=1;
            end
        end
    end
    figure
    imshow(C)
    title('First Image')
    
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
    r=C(:,:,1);
    g=C(:,:,2);
    b=C(:,:,3);
    %% Detections
    r1=215;g1=150;b1=180;
    flag=0;
    while(flag==0)
    im_b= r>r1 & g>g1 & b<b1;
    im_b=imfill(im_b,'holes');
    figure(1),
    subplot(1,2,1)
    imshow(im);
    subplot(1,2,2)
    imshow(im_b);
    choice = questdlg('Satisfied with the detections?','Choose', ...
        'Yes','No','Default');
    cla
    % Handle response
    switch choice
        case 'Yes'
            disp([' Chosen yes'])
            r1=215;g1=150;b1=180;
            flag=1;
        case 'No'
            cla
            figure(2),
            M=imcrop(im);
            temp=mean(mean(M));
            r1=temp(1,1,1);
            b1=temp(1,1,2);
            g1=temp(1,1,3);
    end
    end
    cla
    %% Scoring
    struct1=strel('line',1,30);
    im_b1= r>r1 & g>g1 & b<b1;
    im_b1=imfill(im_b1,'holes');
    figure
    imshow(im_b1)
    title('First Score')
    im_b2= r>r1 & g>g1 & b<b1-20;
    im_b2=imfill(im_b2,'holes');
    figure
    imshow(im_b2)
    title('Second Score')
    im_b3= r>r1 & g>g1 & b<b1-40;
    im_b3=imfill(im_b3,'holes');
    figure
    imshow(im_b3)
    title('Third Score')
    %score=zeros(size(im_b1));
    score=im_b1+im_b2+im_b3;
    figure
    imshow(score)
    title('score')
    [redr,redc]=(find(score==3));
    [or,oc]=(find(score==2));
    [yr,yc]=(find(score==1));
    figure
    subplot(1,2,1)
    imshow(im);
    subplot(1,2,2)
    imshow(im);
    hold on
    plot(yc,yr,'yo');
    hold on
    plot(oc,or,'gx')%'MarkerFaceColor',[ 0.9100 0.4100 0.1700],'MarkerEdgeColor',[ 0.9100 0.4100 0.1700]);
    hold on
    scatter(redc,redr,'r+')%'filled');
    hold on
    %hgexport(gcf, fullfile('plaquescore',filename), hgexport('factorystyle'), 'Format', 'jpeg');
    pause(2)
    %% Process the score matrix to find maximal scores in a circular region around these values
    combined_scoring = score;
    keep_track = score;
    keep_track(:) =0;
    [m,n] = size(score);
    [xgrid, ygrid] = meshgrid(1:n, 1:m);
    
    rad = 10;
   str = strcat('Combined Score within a radius of', ' ', sprintf('%d',rad), ' pixels');
    step_size_x = 10;
    step_size_y = 10;
    [row, col] = find(score > 0);
    min_row = min(row);
    max_row = max(row);
    %centre_pixel_x = 1:step_size_x:n;
    %centre_pixel_y = min_row+rad:step_size_y:max_row-rad;
    for x = 1:step_size_x:n
        for y = min_row+rad:step_size_y:max_row-rad
            mask = ((xgrid-x).^2 + (ygrid-y).^2) <= rad.^2;
            vals = score(mask > 0);
            if any(vals) > 0
                yellow_score = numel(find(vals == 1));
                green_score = numel(find(vals == 2));
                red_score = numel(find(vals == 3));
                ind = score > 0 & mask > 0 & keep_track == 0;
                if yellow_score > green_score && yellow_score > red_score
                    combined_scoring(ind) = 1;
                    keep_track(ind) = 1;
                elseif green_score > yellow_score && green_score > red_score                 
                    combined_scoring(ind) = 2;
                    keep_track(ind) = 1;
                elseif red_score > yellow_score && red_score > green_score
                    combined_scoring(ind) = 3;
                    keep_track(ind) = 1;
                end        
            end
        end
    end
    
    [redr2,redc2]=(find(combined_scoring==3));
    [or2,oc2]=(find(combined_scoring==2));
    [yr2,yc2]=(find(combined_scoring==1));
    figure
    subplot(2,2,1)
    imshow(im);
    title('Original Image')
%     hold on
%     plot(yc,yr,'yo');
%     hold on
%     plot(oc,or,'gx')%'MarkerFaceColor',[ 0.9100 0.4100 0.1700],'MarkerEdgeColor',[ 0.9100 0.4100 0.1700]);
%     hold on
%     scatter(redc,redr,'r+')%'filled');
%     hold on
    subplot(2,2,2)
    imshow(im);
    title('Score = 1 Area')
    hold on
    plot(yc2,yr2,'yo');
    subplot(2,2,3)
    imshow(im);
    title('Score = 2 Area')
    hold on
    plot(oc2,or2,'gx')%'MarkerFaceColor',[ 0.9100 0.4100 0.1700],'MarkerEdgeColor',[ 0.9100 0.4100 0.1700]);
    subplot(2,2,4)
    imshow(im);
    title('Score = 3 Area')
    hold on
    scatter(redc2,redr2,'r+')%'filled');
    hold on       
    %cla
end