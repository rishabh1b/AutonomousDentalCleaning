fileID = fopen('traj.txt','w');
[m,~] = size(upper_teeth_points);
figure 
imshow(im)
hold on
for i = 1 : m
   plot(upper_teeth_points_orig(i,1),upper_teeth_points_orig(i,2),'r*')
   A1=[upper_teeth_points_orig(i,1),upper_teeth_points_orig(i,2), 0.57972,3.1,0.002,1.8415,0,0,0,0,0,0,0,0,0,0,0,0];
   formatSpec = '%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n';
fprintf(formatSpec,A1);
end
fclose(fileID);