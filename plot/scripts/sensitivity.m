M = 100; % Number of macroreplications
num_proc_settings = 4; % Number of processor number settings

% Read in the data from a text file "sensitivity_data.txt"
fileID = fopen('../raw/sensitivity_data.txt','r'); % Currently random data
formatSpec = '%f';
data = fscanf(fileID,formatSpec,[M, num_proc_settings]);
fclose(fileID);

% Extract the obj fn values into column vectors
obj_fn_values_p2 = data(:,1);
obj_fn_values_p4 = data(:,2);
obj_fn_values_p8 = data(:,3);
obj_fn_values_p16 = data(:,4);

% Fix common bin centers across subplots
bin_centers = 0.05:0.1:0.95;

subplot(2,2,1)
hist(obj_fn_values_p2, bin_centers);
title('Histogram of y* for p=2 processors');

subplot(2,2,2);
hist(obj_fn_values_p4, bin_centers);
xlabel('t^* values');
title('Histogram of y* for p=4 processors');

subplot(2,2,3);
hist(obj_fn_values_p8, bin_centers);
xlabel('t^* values');
title('Histogram of y* for p=8 processors');

subplot(2,2,4);
hist(obj_fn_values_p16, bin_centers);
xlabel('t^* values');
title('Histogram of y* for p=16 processors');