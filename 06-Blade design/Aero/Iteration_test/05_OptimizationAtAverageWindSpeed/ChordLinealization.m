function Chord_New = ChordLinealization(Chord_Old, mu)

Len = length(Chord_Old);
Chord_New = zeros(Len,1);

mu_reference = 0.8; 
delta_mu = mu(2)-mu(1);

%Find the chord value at the reference mu point
i_up = find(mu>mu_reference,1);
factor = (mu_reference-mu(i_up-1))/delta_mu;
Chord_reference = Chord_Old(i_up-1)*(1-factor) + Chord_Old(i_up)*factor;

%Find the slope of the chord at this point
m = -(Chord_Old(i_up-1)-Chord_Old(i_up))/delta_mu;

%y = m*x + n -> find n
n = Chord_reference - m * mu_reference;

for i=1:Len
    Chord_New(i) = m*mu(i) + n;
end

end