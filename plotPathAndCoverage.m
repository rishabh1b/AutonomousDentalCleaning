function plotPathAndCoverage(im_mod, pointsForPlotting_u, pointsForPlotting_l, intermediatePoints_u, intermediatePoints_l, laser_rad)
%% Check the path
figure
imshow(im_mod), title('Path by laser pointer')
hold on
% Upper
delta = pointsForPlotting_u(2:end,:) - pointsForPlotting_u(1:end-1,:);
delta = [delta;[0 0]];
quiver(pointsForPlotting_u(:,1),pointsForPlotting_u(:,2),delta(:,1),delta(:,2),0, 'm','MaxHeadSize',0.2, 'LineWidth',2)
% Lower
delta = pointsForPlotting_l(2:end,:) - pointsForPlotting_l(1:end-1,:);
delta = [delta;[0 0]];
quiver(pointsForPlotting_l(:,1),pointsForPlotting_l(:,2),delta(:,1),delta(:,2),0, 'm','MaxHeadSize',0.2, 'LineWidth',2)

% Show Coverage
%pointsForPlotting = [pointsForPlotting, laser_rad .* ones(size(pointsForPlotting,1),1)];
figure 
%im_mod_2 = insertShape(im_mod,'Circle', pointsForPlotting, 'Color', 'green', 'Opacity', 0.6);
imshow(im_mod), title('Coverage by laser')
hold on
s = scatter(intermediatePoints_u(:,1), intermediatePoints_u(:,2),laser_rad*15, 'filled', 'g');
alpha(s, 0.2)
s = scatter(intermediatePoints_l(:,1), intermediatePoints_l(:,2),laser_rad*15, 'filled', 'g');
alpha(s, 0.2)
end