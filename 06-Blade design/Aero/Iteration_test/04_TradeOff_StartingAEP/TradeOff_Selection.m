function Ind = TradeOff_Selection(Blade, Data,Operation)

%AEP Score
AEP_max = max(Blade.AEP);
AEP_min = min(Blade.AEP);
AEP_score = (Blade.AEP-AEP_min)./(AEP_max-AEP_min);

%Starting time score
Ts_max = max(Blade.Ts);
Ts_min = min(Blade.Ts);
Ts_score_old = Ts_min./(Blade.Ts);
Ts_score = (Ts_score_old-min(Ts_score_old))./(1-min(Ts_score_old));

%Thrust Score
for i = 1:length(Operation)
    Thrust_equivalent(i) = trapz(Operation(i).T);
end
Thrust_min = min(Thrust_equivalent);
if(Thrust_min<0)
    Thrust_min = 0;
end
Thrust_score_old = Thrust_min./Thrust_equivalent;
Thrust_score = (Thrust_score_old-min(Thrust_score_old))./(1-min(Thrust_score_old));


Thrust_value = 0.5;

%Total score
Score = AEP_score*Data.a_power_in + Ts_score*(1-Data.a_power_in); % + Thrust_score*Thrust_value;

%Eliminate those options that have stange AoA
% for i = 1:length(Operation)
%     if sum(Operation(i).MeanAoA > 20) >= 1
%         Score(i) = 0;
%     end
% end


%Find the best blade
[~,Ind] = max(Score);
