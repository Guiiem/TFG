function blade_deopt(data)
%function [blade, best_blade, data] = blade_deopt(data)
tic;% record start time
creator;% run creator script to generate initial population
perc_in = 0.05;%  The reporting interval (percentage) 
percent = int8(double(data.no_gen)*perc_in);
fid = fopen(data.fitfile,'w'); %Open output file for fittest blades
disp(' gen.	Max. Cp for blade Min. Start Time for blade')

% calculate fittness of initial population 
blade.fitness = survival(data,blade.chord,blade.twist);

% adjust fitness for poor porforming blades
index = find(blade.fitness(:,2) < data.Cp_min);
if (~isempty(index))
	blade.fitness(index,2) = 0;
end 

% evolve blade population over no_gen generations
for i=1:data.no_gen
	% Find the fittest members of the current population
	[i_max_fit, max_fit, blade,data] = fittest(i,data,blade);
	data.max_fit(i,:) = max_fit;
	if (mod(i,percent)==0)% If timeto report....
		fprintf('%5.1f%% complete,    non-dominated; %5.1f%%   Unfit; %5.1f%%\n',...
            100.0*i/data.no_gen,100*data.no_nondom/data.no_indiv,data.no_unfit/data.no_indiv*100)
		disp(' gen.	Max. Cp for blade Min. Start Time	for blade')
		for j=1:data.no_indiv
			if (~blade.dom(j))
				fprintf(fid,'%12i %12.7f %12.7f \n',j, blade.fitness(j,1), blade.fitness(j,2));
			end
		end	

	end
	if i< data.no_gen	% breed population, except for last gen
		[blade, data] = breed(i, i_max_fit, max_fit,data,blade);
    end
end
fclose(fid); % close file for fittest blades

index = find(blade.dom==0);
% Find the 'best' blade by averaging dominant blades
lindex=length(index);
if lindex >= 1 
    fprintf('\n %3i blades are non-dominated \n \n', length(index))
	best_blade.chord = blade.chord(index,:);
	best_blade.twist = blade.twist(index,:);
elseif isempty(index)
	fprintf(' No blades are non-dominated \n')
end
% prepare to output average of non-dominated blades
best_blade.chord(lindex+1,:)=mean(best_blade.chord);
best_blade.twist(lindex+1,:)=mean(best_blade.twist);
for i = 1:lindex+1
    if (i <= lindex)
        fprintf(' For blade number %5i \n', index(i));
        iprint = 0;
    else
        fprintf(' For the average of the non-dominated blades \n');
        iprint=1;
    end
    best_blade.fitness(i,1)=power_calc_opt(data,best_blade.chord(i,:),best_blade.twist(i,:),iprint);
    J = inertia(best_blade.chord(i,:),best_blade.twist(i,:),data);
    delt = min(data.t_max/50, 0.25);
    best_blade.fitness(i,2)=start_calc_opt(data,best_blade.chord(i,:),best_blade.twist(i,:),J, delt);
    fprintf(' Cp = %5.3f, Starting time = %6.3f seconds \n',best_blade.fitness(i,1), best_blade.fitness(i,2))
    fprintf('\n    radius    chord    twist \n');
    out=[data.rad' best_blade.chord(i,:)' best_blade.twist(i,:)'];
    disp(out)
end
disp(' ');
fprintf('Computation time: %g seconds \n',toc);
end %function blade_deopt.m