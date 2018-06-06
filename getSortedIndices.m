function [true_upper_indices_c_s, true_upper_indices_r_s, true_lower_indices_c_s, true_lower_indices_r_s] = getSortedIndices(plaqueImage)
% Separate in Lower and Upper Points
[m,~] = size(plaqueImage);
[true_indices_r, true_indices_c] = find(plaqueImage);
middle_row = m/2;
logical_upper = true_indices_r < middle_row;
true_upper_indices_r = true_indices_r(logical_upper);
true_upper_indices_c = true_indices_c(logical_upper);
true_lower_indices_r = true_indices_r(~logical_upper);
true_lower_indices_c = true_indices_c(~logical_upper);
% Test
%     figure, imshow(yellow_image), title('Checking Upper Teeth Points')
%     hold on
%     plot(true_upper_indices_c, true_upper_indices_r, 'yo')
%     plot(true_lower_indices_c, true_lower_indices_r, 'rx')

% Sort Upper Points column wise
  [true_upper_indices_c_s, index] = sort(true_upper_indices_c);
  true_upper_indices_r_s = true_upper_indices_r(index);
% Sort Lower Points column wise
  [true_lower_indices_c_s, index] = sort(true_lower_indices_c);
  true_lower_indices_r_s = true_lower_indices_r(index);
end