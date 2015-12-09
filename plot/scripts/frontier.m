% M = # of macroreplications (of the PSCP procedure)
M = 50;

num_proc = [1, 2, 4, 8, 16];
num_settings = length(num_proc);

% Number of assets
num_assets = 200; % num_assets = length(y*)

% Read in the data from a text file "frontier_data.txt"
sizedata = [3+num_assets, M*num_settings];
fileID = fopen('../raw/frontier_data.txt','r');
formatSpec = '%f';
data = fscanf(fileID,formatSpec, sizedata);
fclose(fileID);

data = data';

% Assume the data is sitting in a matrix with columns:
% p (# of processors) | time (sec) | t* | y*

% R = # of realizations (of payouts)
R = 5000;

% Set confidence level for CIs
alpha = 0.05;
z_alpha_over_2 = norminv(1-alpha/2);
t_alpha_over_2 = tinv(1-alpha/2, M-1);

%%
% Suppose we extract all of the times into a matrix "times"
% with M rows and columns of 1|2|4|8|16 processors.
%times = rand(M,num_settings); % fake data set
times = data(:,2);
times = reshape(times,M,num_settings);

% Make plot of wall clock time vs number of processors
% Normality assumption
avg_times = mean(times);
% var_times = var(times);
% lower_CI_times = z_alpha_over_2*(var_times/sqrt(M));
% upper_CI_times = z_alpha_over_2*(var_times/sqrt(M));

figure
plot(log2(num_proc), avg_times);
%errorbar(log2(num_proc), avg_times, lower_CI_times, upper_CI_times);
xlabel('Number of Processors')
ylabel('Wall Clock Time (sec)')
title('Wall Clock Time vs No. of Processors')

V = axis;
V(1:2) = [-0.5, 4.5];
axis(V);

xticks = 0:4;
set(gca, 'XTick', xticks);
xtl = {'1','2','4','8','16'};
set(gca, 'XTickLabel', xtl)

%%
% Suppose we extract all of the obj fn values (t*) into a matrix "objfnvalues"
% with M rows and columns of 1|2|4|8|16 processors.
objfnvalues = data(:,3);
objfnvalues = reshape(objfnvalues,M,num_settings);

% Make plot of obj fn value vs number of processors
% Normality assumption
avg_objfnvalues = mean(objfnvalues);
var_objfnvalues = var(objfnvalues);
lower_CI_objfnvalues = z_alpha_over_2*(var_objfnvalues/sqrt(M));
upper_CI_objfnvalues = z_alpha_over_2*(var_objfnvalues/sqrt(M));

figure
errorbar(log2(num_proc(2:num_settings)), avg_objfnvalues(2:num_settings), lower_CI_objfnvalues(2:num_settings), upper_CI_objfnvalues(2:num_settings));
hold on
plot([1,4], [avg_objfnvalues(1),avg_objfnvalues(1)]);
plot([1,4], [avg_objfnvalues(1)-lower_CI_objfnvalues(1),avg_objfnvalues(1)-lower_CI_objfnvalues(1)],'b:');
plot([1,4], [avg_objfnvalues(1)+upper_CI_objfnvalues(1),avg_objfnvalues(1)+upper_CI_objfnvalues(1)],'b:');
hold off

xlabel('Number of Processors')
ylabel('Obj. Function Value')
title('Obj. Function Value vs No. of Processors')

V = axis;
V(1:2) = [0.5, 4.5];
axis(V);

xticks = 1:4;
set(gca, 'XTick', xticks);
xtl = {'2','4','8','16'};
set(gca, 'XTickLabel', xtl)

%%
% Make plot of obj fn value/wall clock time tradeoff
figure
plot(avg_times, avg_objfnvalues);
text(avg_times, avg_objfnvalues, {'1','2','4','8','16'});
xlabel('Avg Wall Clock Time (sec)');
ylabel('Avg Obj. Function Value');
title('Wall Clock Time / Obj. Function Value Tradeoff');

%%
avg_viol_prob = zeros(1,num_settings);
avg_prob_viol_prob_gr_eps = zeros(1,num_settings);
lower_CI_est_viol_prob = zeros(1,num_settings);
upper_CI_est_viol_prob = zeros(1,num_settings);
lowerbound = zeros(1,num_settings);
upperbound = zeros(1,num_settings);
lower_CI_viol_prob_gr_eps = zeros(1,num_settings);
upper_CI_viol_prob_gr_eps = zeros(1,num_settings);

for i = 1:num_settings

	% Suppose we extract all of the y* and t* into a matrix "solns"
	solns = data((i-1)*M+1:i*M,4:3+num_assets);
    t_values = data((i-1)*M+1:i*M,3);

	% Construct unbiased estimators of E[V(y*)] and Pr(V(y*)>epsilon)
	epsilon = 0.05;

	% This all is just for one number of processors, not all
	% Need to further vectorize
	est_viol_prob = zeros(1,M);
	viol_prob_gr_eps = zeros(1,M);
	for m = 1:M
		ystar = solns(m,1:num_assets);
		%tstar = solns(m,num_assets+1);
		tstar = t_values(m); % for testing

		% Generate a random set of R realizations
		%realizations = rand(num_assets,R);
		realizations = lognrnd(0.094278,sqrt(0.002064),num_assets-1,R);
        realizations = [realizations; 1.05*ones(1,R)];
        
		returns = ystar*realizations;
		est_viol_prob(m) = (1/R)*sum(returns < tstar);
		viol_prob_gr_eps(m) = (est_viol_prob(m) > epsilon);
	end
	avg_viol_prob(i) = mean(est_viol_prob);
	avg_prob_viol_prob_gr_eps(i) = mean(viol_prob_gr_eps);
    x = sum(viol_prob_gr_eps);

	% Construct CIs (Normal approximation of Bernoulli)
	%lower_CI_est_viol_prob(i) = min(z_alpha_over_2*(sqrt(avg_viol_prob(i)*(1-avg_viol_prob(i))/sqrt(M))),avg_viol_prob(i));
	%upper_CI_est_viol_prob(i) = z_alpha_over_2*(sqrt(avg_viol_prob(i)*(1-avg_viol_prob(i))/sqrt(M)));
	%lower_CI_viol_prob_gr_eps(i) = min(z_alpha_over_2*(sqrt(avg_prob_viol_prob_gr_eps(i)*(1-avg_prob_viol_prob_gr_eps(i))/sqrt(M))),avg_prob_viol_prob_gr_eps(i));
	%upper_CI_viol_prob_gr_eps(i) = z_alpha_over_2*(sqrt(avg_prob_viol_prob_gr_eps(i)*(1-avg_prob_viol_prob_gr_eps(i))/sqrt(M)));
    
    % Construct CIs (Normality assumption from CLT)
    lower_CI_est_viol_prob(i) = min(t_alpha_over_2*std(est_viol_prob)/sqrt(M),avg_viol_prob(i));
    upper_CI_est_viol_prob(i) = t_alpha_over_2*std(est_viol_prob)/sqrt(M);
    
    % Construct CIs (Exact Bernoulli CI from Clopper-Pearson
    %%% NEED TO FIX THESE!
    lowerbound(i) = 1/(1+(M-x+1)/(x*finv(alpha/2,2*x,2*(M-x+1))));
    upperbound(i) = 1/(1+(M-x)/((x+1)*finv(1-alpha/2,2*(x+1),2*(M-x))));
    lower_CI_viol_prob_gr_eps(i) = min(avg_prob_viol_prob_gr_eps(i) - lowerbound(i), avg_prob_viol_prob_gr_eps(i));
	upper_CI_viol_prob_gr_eps(i) = upperbound(i) - avg_prob_viol_prob_gr_eps(i);
end

%%
% Make plot of violation probability vs number of processors
figure
errorbar(log2(num_proc(2:num_settings)), avg_viol_prob(2:num_settings), lower_CI_est_viol_prob(2:num_settings), upper_CI_est_viol_prob(2:num_settings));
%errorbar(1, avg_viol_prob, lower_CI_est_viol_prob, upper_CI_est_viol_prob);
xlabel('Number of Processors')
ylabel('E(V(x^*))')
title('E(V(x^*)) vs No. of Processors')

hold on
plot([0.5,4.5], [avg_viol_prob(1),avg_viol_prob(1)],'k');
plot([0.5,4.5], [avg_viol_prob(1)-lower_CI_est_viol_prob(1),avg_viol_prob(1)-lower_CI_est_viol_prob(1)],'k:');
plot([0.5,4.5], [avg_viol_prob(1)+upper_CI_est_viol_prob(1),avg_viol_prob(1)+upper_CI_est_viol_prob(1)],'k:');
hold off

V = axis;
V(1:2) = [0.5, 4.5];
axis(V);

xticks = 1:4;
set(gca, 'XTick', xticks);
xtl = {'2','4','8','16'};
set(gca, 'XTickLabel', xtl)

%%
% Make plot of prob viol prob > epsilon vs number of processors
figure
errorbar(log2(num_proc(2:num_settings)), avg_prob_viol_prob_gr_eps(2:num_settings), lower_CI_viol_prob_gr_eps(2:num_settings), upper_CI_viol_prob_gr_eps(2:num_settings));
%errorbar(1, avg_prob_viol_prob_gr_eps, lower_CI_viol_prob_gr_eps, upper_CI_viol_prob_gr_eps);
xlabel('Number of Processors')
ylabel('Pr(V(x^*) > \epsilon)')
title('Pr(V(x^*) > \epsilon) vs No. of Processors')

hold on
plot([0.5,4.5], [avg_prob_viol_prob_gr_eps(1),avg_prob_viol_prob_gr_eps(1)],'k');
plot([0.5,4.5], [avg_prob_viol_prob_gr_eps(1)-lower_CI_viol_prob_gr_eps(1),avg_prob_viol_prob_gr_eps(1)-lower_CI_viol_prob_gr_eps(1)],'k:');
plot([0.5,4.5], [avg_prob_viol_prob_gr_eps(1)+upper_CI_viol_prob_gr_eps(1),avg_prob_viol_prob_gr_eps(1)+upper_CI_viol_prob_gr_eps(1)],'k:');
hold off

V = axis;
V(1:2) = [0.5, 4.5];
axis(V);

xticks = 1:4;
set(gca, 'XTick', xticks);
xtl = {'2','4','8','16'};
set(gca, 'XTickLabel', xtl)
%%
% Make plot of prob(viol prob > eps)/wall clock time tradeoff
figure
plot(avg_times, avg_prob_viol_prob_gr_eps);
text(avg_times, avg_prob_viol_prob_gr_eps, {'1','2','4','8','16'});
xlabel('Avg Wall Clock Time');
ylabel('Pr(V(x^*) > \epsilon)');
title('Wall Clock Time / Pr(V(x^*) > \epsilon) Tradeoff');
