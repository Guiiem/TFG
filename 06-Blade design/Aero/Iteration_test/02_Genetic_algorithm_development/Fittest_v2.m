function [i_max_fit, max_fit,blade,data] = Fittest_v2(iter, data, blade)

data.no_nondom=0; data.no_unfit=0; data.no_dom=0; 

blade.age = blade.age + 1; 
AEP = blade.AEP;


% store mean and max score for this generation
data.mean_score(iter) = mean(AEP);
data.max_score(iter) = max(AEP);

% Find unfit blades below this generation mean
index = find(AEP < data.mean_score(iter));
if ~isempty(index)
	data.unfit = index;
	blade.dom(index) = deal(true);
 	data.no_unfit = length(index);
    AEP(index)=0;
end

%Find the maximum 
[max_fit,i_max_fit] = max(AEP);

% Find the dominated blades
for i=1:data.no_indiv
	tmp = AEP(i)<AEP;
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

%[0.001*index power(index) time_to_start(index)] ??
% if (iter == data.no_gen)
%     figure(1)
%     plot(time_to_start(index),AEP(index),'x')
% end

% Find number of non-dominated blades
if ~isempty(index)
	data.no_nondom= length(index);
	data.non_dom = index;
end

% disp(sprintf('%4i %8.3f      %4i    %8.3f          %4i', iter, max_fit, i_max_fit))
blade.fitness = AEP;

end