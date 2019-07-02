function [i_max_fit, max_fit,blade,data] = fittest(iter, data, blade)

data.no_nondom=0; data.no_unfit=0; data.no_dom=0; 
max_time = 0.0; min_power = 0.593;

blade.age = blade.age + 1;
power = blade.fitness(:,1);
time_to_start = blade.fitness(:,2);

min_power = min(power);
if ( min_power < data.Cp_min);   min_power = data.Cp_min; end
max_time = max(time_to_start);

% store mean and max score for this generation
data.mean_score(iter) = mean(power);
data.max_score(iter) = max(power);
% Find unfit blades of low power and/or high starting times
index = find(power< data.Cp_min | time_to_start > data.t_max);
if ~isempty(index)
	data.unfit = index;
	blade.dom(index) = deal(true);
 	data.no_unfit = length(index);
    power(index)=data.Cp_min;
    time_to_start(index)=data.t_max;   
end
[max_fit(1),i_max_fit(1)] = max(power);
[max_fit(2),i_max_fit(2)] = min(time_to_start);
% Find the dominated blades
for i=1:data.no_indiv
	tmp = power(i)<power& time_to_start(i)>time_to_start;
	blade.dom(i)=any(tmp);
end
% Find number of dominated blades
index = find(blade.dom);
if ~isempty(index)
	data.no_dom= length(index);
	data.dominated = index;
end

%find non-dominated blades
index = find(~blade.dom);
%[0.001*index power(index) time_to_start(index)]
if (iter == data.no_gen)
    figure(1)
    plot(time_to_start(index),power(index),'x')
end
% Find number of non-dominated blades
if ~isempty(index)
	data.no_nondom= length(index);
	data.non_dom = index;
end

disp(sprintf('%4i %8.3f      %4i    %8.3f          %4i', iter, max_fit(1), i_max_fit(1), max_fit(2), i_max_fit(2)))
blade.fitness(:,1)=power;
blade.fitness(:,2)=time_to_start;

end