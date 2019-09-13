function [blade, data] = Breed_v2(iter, i_max_fit, max_fit, data, blade, plots)
%  Uses the best individuals to breed another population

num_indiv = data.no_indiv; %extract from data structure
num_genes = data.no_genes;

%Calculate fitness score
score = blade.fitness/max_fit;
	
%Initialise Populations
old_chord = blade.Chord;
old_twist = blade.Twist;
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


% LINE TO BE MODIFIED
%trial_fitness = survival(data,trial_chord,trial_twist); %%D'aquí hauria de sortir AEP
for i=1:num_indiv
    [trial_fitness(i), ~, ~] = BigSolver(data, trial_chord(i,:), trial_twist(i,:), plots);
end

%trial_score = data.a_power_in/max_fit(1)*trial_fitness(:,1) + ...
             %   (1-data.a_power_in)*max_fit(2)./trial_fitness(:,2);
trial_score = trial_fitness/max_fit;

% Substitute trial blade for unfit or old blades
for i=1:num_indiv
    if (trial_score(i) > score(i) || blade.age(i) == data.max_age)
        blade.Chord(i,:) = trial_chord(i,:);
        blade.Twist(i,:) = trial_twist(i,:);
        blade.fitness(1,i) = trial_fitness(1,i);
        blade.dom(i) = deal(true);
        blade.age(i) = 1;
    end
end
%index = find(trial_score > score | blade.age.' == data.max_age );

%if (~isempty(index))
% 	blade.Chord(index,:) = trial_chord(index,:);
% 	blade.Twist(index,:) = trial_twist(index,:);
% 	blade.fitness(index,:) = trial_fitness(index,:);
% 	blade.dom(index) = deal(true);
% 	blade.age(index) = 1;
end