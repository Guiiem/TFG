function Blade = Creator (data)

Blade.chord = zeros(data.no_indiv,data.no_genes);
Blade.twist = zeros(data.no_indiv,data.no_genes);

Blade.fitness = zeros(data.no_indiv,2);
Blade.dom = repmat(true,data.no_indiv,1);
Blade.age = ones(data.no_indiv,1);

delc = data.max_chord - data.min_chord;
delt = data.max_twist - data.min_twist;

for i = 1:data.no_indiv    % Generate the first population
    rnd1 = rand;% The blade chord and twist are random, but are
    rnd2 = rand;% constrained to lie between max and min values
    for j = 1:data.no_genes
        Blade.Chord(i,j) = rnd1*delc + data.min_chord;
        Blade.Twist(i,j) = rnd2*delt + data.min_twist;
    end
end



end

