function [blade, data] = breed(iter, i_max_fit, max_fit, data, blade)
%  Uses the best individuals to breed another population

num_indiv = data.no_indiv; %extract from data structure
num_genes = data.no_genes;

%calculate fitness score
score = data.a_power_in/max_fit(1)*blade.fitness(:,1) + ...
(1-data.a_power_in)*max_fit(2)./blade.fitness(:,2);

score(i_max_fit(1)) = 1.0;
score(i_max_fit(2)) = 1.0;
	
%Initialise Populations
old_chord = blade.chord;
old_twist = blade.twist;
init=zeros(num_indiv,num_genes);
basis_chord = init;basis_twist = init; u_chord = init; u_twist = init;
l_chord = init;l_twist = init;int_mask = init; inv_mask = init;
rotating_index = (0:1:num_indiv-1);
init = zeros(num_indiv,1);
rotated =init; basis_index = init; u_index = init; l_index = init;
%calculate mask to for gene swapping
int_mask = rand(num_indiv,num_genes) < data.crossover;		
basis_index = randperm(num_indiv);
rotated = rem(rotating_index+1,num_indiv);
u_index = basis_index(rotated+1);
l_index = u_index(rotated+1);
basis_chord = old_chord(basis_index,:);
basis_twist = old_twist(basis_index,:);
	
u_chord = old_chord(u_index,:); u_twist = old_twist(u_index,:);	
l_chord = old_chord(l_index,:);l_twist = old_twist(l_index,:);
% Generate trial blades	
trial_chord = basis_chord + 0.8*(u_chord - l_chord).*int_mask;		
trial_twist = basis_twist + 0.8*(u_twist - l_twist).*int_mask;	
% Constrain twist and chord of new blades	
trial_chord = min(trial_chord, data.max_chord);
trial_chord = max(trial_chord, data.min_chord);
trial_twist = min(trial_twist, data.max_twist);
trial_twist = max(trial_twist, data.min_twist);
trial_fitness = survival(data,trial_chord,trial_twist);

trial_score = data.a_power_in/max_fit(1)*trial_fitness(:,1) + ...
                (1-data.a_power_in)*max_fit(2)./trial_fitness(:,2);
% Substitute trial blade for unfit or old blades
index = find(trial_score > score | blade.age == data.max_age );
if (~isempty(index))
	blade.chord(index,:) = trial_chord(index,:);
	blade.twist(index,:) = trial_twist(index,:);
	blade.fitness(index,:) = trial_fitness(index,:);
	blade.dom(index) = deal(true);
	blade.age(index) = 1;
end