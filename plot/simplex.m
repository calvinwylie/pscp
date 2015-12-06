% M = # of macroreplications (of the simplex method)
M = 20;

% Read in the data from a text file "simplex_data.txt"
sizedata = [3, M*5];
fileID = fopen('simplex_data.txt','r');
formatSpec = '%f %f %f';
data = fscanf(fileID,formatSpec,sizedata);
fclose(fileID);

% Assume the data is sitting in a matrix with columns:
%  # constraints | simplex time (sec) | dual simplex times (sec)
num_constr = [250, 500, 1000, 2500, 5000]; % # constraints
num_settings = length(num_constr)

avg_simp_times = zeros(num_settings,1);
avg_dual_simp_times = zeros(num_settings,1);

for i = 1:num_settings;
	avg_simp_times(i) = mean(data(2,(i-1)*M+1:i*M));
	avg_dual_simp_times(i) = mean(data(3,(i-1)*M+1:i*M));
end

% Make of a plot of the performance of the two methods (as # constraints increases)
plot(num_constr, avg_simp_times, num_constr, avg_dual_simp_times);
legend('Simplex','Dual Simplex');
xlabel('Number of Constraints')
ylabel('Wall Clock Time (sec)')
title('Time Scaling of Simplex and Dual Simplex Methods')

% Can always try some confidence intervals too. May not be normal.