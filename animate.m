figure(3)
imshow(C)
hold on

for k = 1 : sz_1
    temp_sz = size(selected_cols_upper_sorted_cleaned{k},2);
    for t = 1 : temp_sz
        plot(selected_cols_upper_sorted_cleaned{k}(t),selected_rows_upper_sorted_cleaned{k}(t), 'ro');
        pause(0.2)
    end
end

for k = 1 : sz_2
    temp_sz = size(selected_cols_lower_sorted_cleaned{k},2);
    for t = 1 : temp_sz
        plot(selected_cols_lower_sorted_cleaned{k}(t),selected_rows_lower_sorted_cleaned{k}(t), 'ro');
        pause(0.2)
    end
end