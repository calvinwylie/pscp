% Read in the data from a text file "frontier_data.txt"
fileID = fopen('frontier_data.txt','r');
formatSpec = '%f';
data = fscanf(fileID,formatSpec);
fclose(fileID);

% Assume the data is sitting in a matrix with columns:
%  p (# of processors) | time (sec) | t* | y*
num_proc = [1, 2, 4, 8, 16];

% M = # of macroreplications (of the PSCP procedure)
% R = # of realizations (of payouts)
M = 100;
R = 1000;

% Set confidence level for CIs
alpha = 0.05;
z_alpha_over_2 = norminv(1-alpha/2);

% Suppose we extract all of the times into a matrix "times"
% with M rows and columns of 1|2|4|8|16 processors.
times = rand(M,length(num_proc)); % fake data set

% Make plot of wall clock time vs number of processors
% Normality assumption
avg_times = mean(times);
var_times = var(times);
lower_CI_times = z_alpha_over_2*(var_times/sqrt(M));
upper_CI_times = z_alpha_over_2*(var_times/sqrt(M));
errorbar(log(num_proc), avg_times, lower_CI_times, upper_CI_times);
xlabel('log_2(Number of Processors)')
ylabel('Wall Clock Time (sec)')
title('Wall Clock Time vs No. of Processors')

pause;

% Suppose we extract all of the obj fn values (t*) into a matrix "objfnvalues"
% with M rows and columns of 1|2|4|8|16 processors.
objfnvalues = rand(M,length(num_proc)); % fake data set

% Make plot of obj fn value vs number of processors
% Normality assumption
avg_objfnvalues = mean(objfnvalues);
var_objfnvalues = var(objfnvalues);
lower_CI_objfnvalues = z_alpha_over_2*(var_objfnvalues/sqrt(M));
upper_CI_objfnvalues = z_alpha_over_2*(var_objfnvalues/sqrt(M));
errorbar(log(num_proc), avg_objfnvalues, lower_CI_objfnvalues, upper_CI_objfnvalues);
xlabel('log_2(Number of Processors)')
ylabel('Obj. Function Value')
title('Obj. Function Value vs No. of Processors')

pause;

% Make plot of obj fn value/wall clock time tradeoff
plot(avg_times, avg_objfnvalues);
xlabel('Avg Wall Clock Time');
ylabel('Avg Obj. Function Value');
title('Wall Clock Time / Obj. Function Value Tradeoff');
pause;

% Number of assets
num_assets = 200; % num_assets = length(y*)

% Suppose we extract all of the y* and t* into a matrix "solns"
solns = rand(M,201);

% Construct unbiased estimators of E[V(y*)] and Pr(V(y*)>epsilon)
epsilon = 0.25; % Will be much smaller

% This all is just for one number of processors, not all
% Need to further vectorize
est_viol_prob = zeros(1,M);
viol_prob_gr_eps = zeros(1,M);
for m = 1:M
	ystar = solns(m,1:num_assets);
	%tstar = solns(m,num_assets+1);
	tstar = 50; % for testing

	% Generate a random set of R realizations
	realizations = rand(num_assets,R);
	
	returns = ystar*realizations;
	est_viol_prob(m) = (1/R)*sum(returns > tstar);
	viol_prob_gr_eps(m) = (est_viol_prob(m) > epsilon);
end
avg_viol_prob = mean(est_viol_prob);
avg_prob_viol_prob_gr_eps = mean(viol_prob_gr_eps);

% Construct CIs
lower_CI_est_viol_prob = z_alpha_over_2*(sqrt(avg_viol_prob*(1-avg_viol_prob)/sqrt(M)));
upper_CI_est_viol_prob = z_alpha_over_2*(sqrt(avg_viol_prob*(1-avg_viol_prob)/sqrt(M)));
lower_CI_viol_prob_gr_eps = z_alpha_over_2*(sqrt(avg_prob_viol_prob_gr_eps*(1-avg_prob_viol_prob_gr_eps)/sqrt(M)));
upper_CI_viol_prob_gr_eps = z_alpha_over_2*(sqrt(avg_prob_viol_prob_gr_eps*(1-avg_prob_viol_prob_gr_eps)/sqrt(M)));

% Make plot of violation probability vs number of processors
%errorbar(log(num_proc), avg_viol_prob, lower_CI_est_viol_prob, upper_CI_est_viol_prob);
errorbar(1, avg_viol_prob, lower_CI_est_viol_prob, upper_CI_est_viol_prob);
xlabel('log_2(Number of Processors)')
ylabel('Exp. Violation Prob.')
title('Exp. Violation Prob. vs No. of Processors')

pause;

% Make plot of prob viol prob > epsilon vs number of processors
%errorbar(log(num_proc), avg_prob_viol_prob_gr_eps, lower_CI_viol_prob_gr_eps, upper_CI_viol_prob_gr_eps);
errorbar(1, avg_prob_viol_prob_gr_eps, lower_CI_viol_prob_gr_eps, upper_CI_viol_prob_gr_eps);
xlabel('log_2(Number of Processors)')
ylabel('Prob. Violation Prob. > epsilon')
title('Prob. Violation Prob. > epsilon vs No. of Processors')