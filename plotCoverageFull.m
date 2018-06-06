function plotCoverageFull(im, intermediatePoints_u_y, intermediatePoints_l_y, intermediatePoints_u_b, intermediatePoints_l_b, intermediatePoints_u_r, intermediatePoints_l_r, gumLinepointsLower, gumLinePointsUpper, laser_rad)
global count;

% Show Coverage
%pointsForPlotting = [pointsForPlotting, laser_rad .* ones(size(pointsForPlotting,1),1)];
figure4 = figure(4); 
%im_mod_2 = insertShape(im_mod,'Circle', pointsForPlotting, 'Color', 'green', 'Opacity', 0.6);
imshow(im), title('Coverage by laser')
hold on
% GumLine
s = scatter(gumLinepointsLower(:,1), gumLinepointsLower(:,2),laser_rad*15, 'filled', 'c');
alpha(s, 0.2)
s = scatter(gumLinePointsUpper(:,1), gumLinePointsUpper(:,2),laser_rad*15, 'filled', 'c');
alpha(s, 0.2)
% Yellow
if ~isempty(intermediatePoints_u_y)
    s2 = scatter(intermediatePoints_u_y(:,1), intermediatePoints_u_y(:,2),laser_rad*15, 'filled', 'y');
    alpha(s2, 0.2)
end
if ~isempty(intermediatePoints_l_y)
    s2 = scatter(intermediatePoints_l_y(:,1), intermediatePoints_l_y(:,2),laser_rad*15, 'filled', 'y');
    alpha(s2, 0.2)
end
% Green
if ~isempty(intermediatePoints_u_b)
    s3 = scatter(intermediatePoints_u_b(:,1), intermediatePoints_u_b(:,2),laser_rad*15, 'filled', 'g');
    alpha(s3, 0.3)
end
if ~isempty(intermediatePoints_l_b)
    s3 = scatter(intermediatePoints_l_b(:,1), intermediatePoints_l_b(:,2),laser_rad*15, 'filled', 'g');
    alpha(s3, 0.3)
end
% Red
if ~isempty(intermediatePoints_u_r)
    s4 = scatter(intermediatePoints_u_r(:,1), intermediatePoints_u_r(:,2),laser_rad*15, 'filled', 'r');
    alpha(s4, 0.4)
end
if ~isempty(intermediatePoints_l_r)
    s4 = scatter(intermediatePoints_l_r(:,1), intermediatePoints_l_r(:,2),laser_rad*15, 'filled', 'r');
    alpha(s4, 0.4)
end
if exist('s4', 'var') == 1
    lgd = legend([s, s2, s3, s4], 'GumLine', 'Score1', 'Score2', 'Score3');
else
    lgd = legend([s, s2, s3], 'GumLine', 'Score1', 'Score2');
end
lgd.Location = 'southeast';
lgd.FontSize = 15;
filename = sprintf('%d.jpg', count);
filename = strcat('coverageImages/', filename);
saveas(figure4,filename)
end