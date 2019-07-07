function [Chord,Twist] = ChordTwist(Data);

N = Data.no_genes; %Number of blade elements 
Cmax = Data.max_chord; %Maximum chord
Cmin = Data.min_chord; %Minimum chord
Tmax = Data.max_twist; %Maximum twist
Tmin = Data.min_twist; %Minimum twist

%Chord = fliplr(Cmin:((Cmax-Cmin)/(N-1)):Cmax);
%Twist = fliplr(Tmin:((Tmax-Tmin)/(N-1)):Tmax);

delc = Data.max_chord - Data.min_chord;
delt = Data.max_twist - Data.min_twist;

rnd1 = rand;% The blade chord and twist are random, but are 
rnd2 = rand;% constrained to lie between max and min values 

for i=1:N
    Chord(i) = rnd1*delc + Data.min_chord;
    Twist(i) = rnd2*delt + Data.min_twist;
end

end

