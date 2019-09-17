function Lambda_selected_idx = FindLambdaGeometry(Data, Power, Omega)

len = length(Omega);

%First, calculate the power of the generator at each Omega value of the
%rotor

Torque_Gen = Omega.*Data.gen_torque(2)/Data.gen_omega(2);

Torque_Rotor = Power./(Omega*2*pi/60);

%Localize the minimum difference
for i = 1:len
    Diff(i) = abs(Torque_Gen(i)-Torque_Rotor(i));
end

[~, Lambda_selected_idx] = min(Diff(:));


end

