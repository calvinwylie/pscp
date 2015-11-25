% Read in the data from a text file "sensitivity_data.txt"
fileID = fopen('sensitivity_data.txt','r');
formatSpe = '%f';
data = fscanf(fileID,formatSpec);

% Extract the obj fn values into column vectors
obj_fn_values_p2 = data(:,1);
obj_fn_values_p4 = data(:,2);
obj_fn_values_p8 = data(:,3);
obj_fn_values_p16 = data(:,4);

% Fix common bin centers across subplots
%bin_centers = [];

subplot(2,2,1)
hist(obj_fn_values_p2);
xlabel('t^* values');
title('Histogram of y* for p=2 processors');

suplot(2,2,2);
hist(obj_fn_values_p4);
xlabel('t^* values');
title('Histogram of y* for p=4 processors');

subplot(2,2,3);
hist(obj_fn_values_p8);
xlabel('t^* values');
title('Histogram of y* for p=8 processors');

subplot(2,2,4);
hist(obj_fn_values_p16);
xlabel('t^* values');
title('Histogram of y* for p=16 processors');