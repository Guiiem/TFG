clear variables

lambda = linspace(0,10);
U = 10; %ms
R = 0.6; %m
omega1 = lambda.*U/R;
omega2 = omega1*60/(2*pi);
d = 0.3;
m = 0.4; %kg
F = m*d*omega1.^2;

figure()
yyaxis left
plot(lambda,omega2)
xlabel('Tip speed ratio [-]')
ylabel('Rotor speed [rpm]')

yyaxis right
plot(lambda,F)
ylabel('Centrifugal force [N]')
