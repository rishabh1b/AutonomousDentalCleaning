show_image = true;
global count
cleaning_times = [];
for i = 48:1:48
    count = i;
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
    if show_image
        figure
        imshow(C)
        title('First Image')
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
%     choice = questdlg('Satisfied with the detections?','Choose', ...
%         'Yes','No','Default');
    choice = 'Yes';
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
    if show_image
        figure
        imshow(im_b1)
        title('First Score')
    end
    im_b2= r>r1 & g>g1 & b<b1-20;
    im_b2=imfill(im_b2,'holes');
    if show_image
        figure
        imshow(im_b2)
        title('Second Score')
    end
    im_b3= r>r1 & g>g1 & b<b1-40;
    im_b3=imfill(im_b3,'holes');
    if show_image
        figure
        imshow(im_b3)
        title('Third Score')
    end
    %score=zeros(size(im_b1));
    score=im_b1+im_b2+im_b3;
    if show_image
        figure
        imshow(score)
        title('score')
    end
    [redr,redc]=(find(score==3));
    [or,oc]=(find(score==2));
    [yr,yc]=(find(score==1));
    if show_image
        figure
        imshow(im)
        hold on
        plot(yc,yr,'yo')
        plot(oc,or,'gx')
        scatter(redc,redr,'r+')
        lgd = legend('Score 1', 'Score 2', 'Score 3');
        lgd.FontSize = 20;
%         subplot(1,2,1)
%         imshow(im);
%         subplot(1,2,2)
%         imshow(im);
%         hold on
%         plot(yc,yr,'yo');
%         hold on
%         plot(oc,or,'gx')%'MarkerFaceColor',[ 0.9100 0.4100 0.1700],'MarkerEdgeColor',[ 0.9100 0.4100 0.1700]);
%         hold on
%         scatter(redc,redr,'r+')%'filled');
%         hold on
%         %hgexport(gcf, fullfile('plaquescore',filename), hgexport('factorystyle'), 'Format', 'jpeg');
%         pause(2)
    end
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
    if show_image
        figure
        subplot(2,2,1)
        imshow(im);
        title('Original Image')
        subplot(2,2,2)
        imshow(im);
        title('Score = 1 Area')
        hold on
        plot(yc2,yr2,'yo');
        subplot(2,2,3)
        imshow(im);
        title('Score = 2 Area')
        hold on
        plot(oc2,or2,'bx')%'MarkerFaceColor',[ 0.9100 0.4100 0.1700],'MarkerEdgeColor',[ 0.9100 0.4100 0.1700]);
        subplot(2,2,4)
        imshow(im);
        title('Score = 3 Area')
        hold on
        scatter(redc2,redr2,'r+')%'filled');
        hold on     
    end
    %cla
    %% Create some binary images
    red_image = false(m,n);
    green_image = false(m,n);
    yellow_image = false(m,n);
    num_rows = size(redc2,1);
    for j = 1 : num_rows
        red_image(redr2(j), redc2(j)) = true;
    end
    num_rows = size(or2,1);
    for j = 1 : num_rows
        green_image(or2(j), oc2(j)) = true;
    end
    num_rows = size(yr2,1);
    for j = 1 : num_rows
        yellow_image(yr2(j), yc2(j)) = true;
    end
    % Remove small noises in the middle
    red_image = bwareaopen(red_image, 20);
    green_image = bwareaopen(green_image, 20);
    yellow_image = bwareaopen(yellow_image, 20);
    if show_image
        figure
        imshow(red_image), title('Red Plaque Regions')
        figure
        imshow(green_image), title('Green Plaque Regions')
        figure
        imshow(yellow_image), title('Yellow Plaque Regions')
    end
    %% Get the Indices in sorted order for different plaque images
    [true_upper_indices_c_s_y, true_upper_indices_r_s_y, true_lower_indices_c_s_y, true_lower_indices_r_s_y] = getSortedIndices(yellow_image);
    [true_upper_indices_c_s_b, true_upper_indices_r_s_b, true_lower_indices_c_s_b, true_lower_indices_r_s_b] = getSortedIndices(green_image);
    [true_upper_indices_c_s_r, true_upper_indices_r_s_r, true_lower_indices_c_s_r, true_lower_indices_r_s_r] = getSortedIndices(red_image);
%     % Separate in Lower and Upper Points
%     [true_indices_r, true_indices_c] = find(yellow_image);
%     middle_row = m/2;
%     logical_upper = true_indices_r < middle_row;
%     true_upper_indices_r = true_indices_r(logical_upper);
%     true_upper_indices_c = true_indices_c(logical_upper);
%     true_lower_indices_r = true_indices_r(~logical_upper);
%     true_lower_indices_c = true_indices_c(~logical_upper);
%     % Test
% %     figure, imshow(yellow_image), title('Checking Upper Teeth Points')
% %     hold on
% %     plot(true_upper_indices_c, true_upper_indices_r, 'yo')
% %     plot(true_lower_indices_c, true_lower_indices_r, 'rx')
% 
%     % Sort Upper Points column wise
%       [true_upper_indices_c_s, index] = sort(true_upper_indices_c);
%       true_upper_indices_r_s = true_upper_indices_r(index);
%     % Sort Lower Points column wise
%       [true_lower_indices_c_s, index] = sort(true_lower_indices_c);
%       true_lower_indices_r_s = true_lower_indices_r(index);
    %% Blend (Red, Yellow and Blue) coloured plaque on top of the teeth
    % Add Transparent Screen of yellow on top of the image
    im_mod_y = im;
    im_r = im(:,:,1);
    im_g = im(:,:,2);
    im_b = im(:,:,3);
    im_g_new = uint8(zeros(size(im,1), size(im,2)));
    im_r_new = im_g_new;
    im_r_new(yellow_image == 1) = 154;
    im_g_new(yellow_image == 1) = 100;
    im_b(yellow_image == 1) = 0;
    im_g_2 = im_g;
    im_r_2 = im_r;
    % Blend the images to get a transparent effect
    im_g_2(yellow_image == 1) = 0.8 * im_g(yellow_image == 1) + 0.2 * im_g_new(yellow_image == 1);
    im_r_2(yellow_image == 1) = 0.8 * im_r(yellow_image == 1) + 0.4 * im_r_new(yellow_image == 1);
    im_mod_y(:,:,1) = im_r_2;
    im_mod_y(:,:,2) = im_g_2;
    im_mod_y(:,:,3) = im_b;
    if show_image
        figure
        imshow(im_mod_y), title('Yellow Film Overlayed')
    end
    
     % Add Transparent Screen of Blue on top of the image
    im_mod_b = im;
    im_r = im(:,:,1);
    im_g = im(:,:,2);
    im_b = im(:,:,3);
    im_b_new = uint8(zeros(size(im,1), size(im,2)));
    im_r(green_image == 1) = 0;
    im_g(green_image == 1) = 0;
    im_b_new(green_image == 1) = 127;
    im_b_2 = im_b;
    % Blend the images to get a transparent effect
    im_b_2(green_image == 1) = 0.8 * im_b(green_image == 1) + 0.02 * im_b_new(green_image == 1);
    im_mod_b(:,:,1) = im_r;
    im_mod_b(:,:,2) = im_g;
    im_mod_b(:,:,3) = im_b_2;
    if show_image
        figure
        imshow(im_mod_b), title('Blue Film Overlayed')
    end
    
    % Add Transparent Screen of Red on top of the image
    im_mod_r = im;
    im_r = im(:,:,1);
    im_g = im(:,:,2);
    im_b = im(:,:,3);
    im_r_new = uint8(zeros(size(im,1), size(im,2)));
    im_r_new(red_image == 1) = 255;
    im_g(red_image == 1) = 0;
    im_b(red_image == 1) = 0;
    im_r_2 = im_r;
    % Blend the images to get a transparent effect
    im_r_2(red_image == 1) = 0.8 * im_r(red_image == 1) + 0.02 * im_r_new(red_image == 1);
    im_mod_r(:,:,1) = im_r_2;
    im_mod_r(:,:,2) = im_g;
    im_mod_r(:,:,3) = im_b;
    if show_image
        figure
        imshow(im_mod_r), title('Red Film Overlayed')
    end
    %% Call the algorithm
    laser_rad = 5;
    intermediatePoints_u_y = [];
    intermediatePoints_l_y = [];
    intermediatePoints_u_b = [];
    intermediatePoints_l_b = [];
    intermediatePoints_u_r = [];
    intermediatePoints_l_r = [];
    % Yellow Teeth
    if ~isempty(true_upper_indices_c_s_y) 
        [pointsForPlotting_u_y, intermediatePoints_u_y] = coverageMapping(true_upper_indices_c_s_y, true_upper_indices_r_s_y, 2, laser_rad);
    end
    if ~isempty(true_lower_indices_c_s_y)
        [pointsForPlotting_l_y, intermediatePoints_l_y] = coverageMapping(true_lower_indices_c_s_y, true_lower_indices_r_s_y, 1, laser_rad);
    end
    %plotPathAndCoverage(im_mod_y, pointsForPlotting_u, pointsForPlotting_l, intermediatePoints_u, intermediatePoints_l, laser_rad)
    % Green/Blue Regions
     if ~isempty(true_upper_indices_c_s_b)
        [pointsForPlotting_u_b, intermediatePoints_u_b] = coverageMapping(true_upper_indices_c_s_b, true_upper_indices_r_s_b, 2, laser_rad);
     end
     if ~isempty(true_lower_indices_c_s_b)
        [pointsForPlotting_l_b, intermediatePoints_l_b] = coverageMapping(true_lower_indices_c_s_b, true_lower_indices_r_s_b, 1, laser_rad);
     end
    %plotPathAndCoverage(im_mod_b, pointsForPlotting_u, pointsForPlotting_l, intermediatePoints_u, intermediatePoints_l, laser_rad)
    % Red Regions
     if ~isempty(true_upper_indices_c_s_r)
         [pointsForPlotting_u_r, intermediatePoints_u_r] = coverageMapping(true_upper_indices_c_s_r, true_upper_indices_r_s_r, 2, laser_rad);
     end
     if ~isempty(true_lower_indices_c_s_r)
        [pointsForPlotting_l_r, intermediatePoints_l_r] = coverageMapping(true_lower_indices_c_s_r, true_lower_indices_r_s_r, 1, laser_rad);
     end
    %plotPathAndCoverage(im_mod_r, pointsForPlotting_u, pointsForPlotting_l, intermediatePoints_u, intermediatePoints_l, laser_rad)
    
    % Get the Points on Gum Line
    getToothProfile2;
    %% Plot Entire Coverage
    plotCoverageFull(im, intermediatePoints_u_y, intermediatePoints_l_y, intermediatePoints_u_b, intermediatePoints_l_b, intermediatePoints_u_r, intermediatePoints_l_r, gumLinepointsLower, gumLinePointsUpper, laser_rad)
    %% Evaluate the coverage time
    timeForCleaning = size(gumLinePointsUpper,1) / 7 + size(gumLinepointsLower,1) / 7 +...
                       size(intermediatePoints_u_y,1) / 7 + size(intermediatePoints_l_y,1) / 7 + ...
                       size(intermediatePoints_u_b,1) / 5 + size(intermediatePoints_l_b,1) / 5 + ...
                       size(intermediatePoints_u_r,1) / 3 + size(intermediatePoints_l_r,1) / 3;
   fprintf('The time for cleaning for Image %d is %f min \n', i, timeForCleaning/60);
   cleaning_times = [cleaning_times; timeForCleaning/60];
end