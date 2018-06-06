% Init Algorithm for map creation
%TODO - Think about shift in pixels correlation to metric
%TODO - Incorporate Time metrices in the code
% Intiate everything in unvisited array and dontconsider array
num_points = size(true_upper_indices_c_s,1);
unvisited_c = true_upper_indices_c_s;
unvisited_r = true_upper_indices_r_s;
dontconsider = false(size(true_upper_indices_c_s,1),1);
% Initate a checkpost points array
intermediatePoints = [];
pointsForPlotting = [];
% Initiate shifting array
shift_c = 5;
laser_rad = 10;
shift_r = 5;
shifting_array = [shift_r, -shift_r];
shifting_ind = 2;
% Initiate current postion in the unvisited set
curr_ind = 1;
% Take the first pixel value in the array and assign it to centre
cx = unvisited_c(curr_ind);
cy = unvisited_r(curr_ind);
pointsForPlotting = [pointsForPlotting;[cx cy]];
while any(~dontconsider) || dontconsider(end) == true
    currSubset_c = unvisited_c(curr_ind:end,1);
    currSubset_r = unvisited_r(curr_ind:end,1);
    pointcheck = ((currSubset_c-cx).^2 + (currSubset_r-cy).^2) <= laser_rad.^2;
    pointcheck = pointcheck & (~dontconsider);
    while any(pointcheck)
        dontconsider(pointcheck > 0) = true;
        %curr_ind = curr_ind + 1;
        intermediatePoints = [intermediatePoints;[cx cy]];
        cy = cy  + shifting_array(shifting_ind);
        currSubset_c = unvisited_c(curr_ind:end,1);
        currSubset_r = unvisited_r(curr_ind:end,1);
        pointcheck = ((currSubset_c-cx).^2 + (currSubset_r-cy).^2) <= laser_rad.^2;
    end
% Make a Right and check whether any 'new' pixel value can be found
    cx = cx + shift_c;
    pointcheck = ((currSubset_c-cx).^2 + (currSubset_r-cy).^2) <= laser_rad.^2;
    %pointcheck = pointcheck & (~dontconsider);
    if ~any(pointcheck) % Check whether we have points in this direction 
        % Otherwise go for the next unvisited point in the array
        ind = find(dontconsider);
        if isempty(ind) || ind(end) == num_points
            break;
        end
        cx = unvisited_c(ind(end) + 1);
        cy = unvisited_r(ind(end) + 1);
    end

    % Change the shifting style
    pointsForPlotting = [pointsForPlotting;[cx cy]];
    if shifting_ind == 2
        shifting_ind = 1;
    else
        shifting_ind = 2;
    end
end
%% Check the path
figure
imshow(im_mod), title('Path by laser pointer')
hold on
delta = pointsForPlotting(2:end,:) - pointsForPlotting(1:end-1,:);
delta = [delta;[0 0]];
quiver(pointsForPlotting(:,1),pointsForPlotting(:,2),delta(:,1),delta(:,2),0, 'MaxHeadSize',0.2)
%plot(intermediatePoints(:,1), intermediatePoints(:,2), 'b*', 'MarkerSize', 2)

% Show Coverage
%pointsForPlotting = [pointsForPlotting, laser_rad .* ones(size(pointsForPlotting,1),1)];
figure 
%im_mod_2 = insertShape(im_mod,'Circle', pointsForPlotting, 'Color', 'green', 'Opacity', 0.6);
imshow(im_mod), title('Coverage by laser')
hold on
plot(intermediatePoints(:,1), intermediatePoints(:,2), 'go', 'MarkerSize', 15, 'MarkerFaceColor', 'green')


%plot(true_upper_indices_c_s, true_upper_indices_r_s, 'b*')