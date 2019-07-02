function [Qeq Omegaeq] = Intersect(Q,Omega,QGen,OmegaGen,U);

Qeq=NaN;
Omegaeq=NaN;

%First, compute the Qgen at the discrete Omega points
for i=1:length(Omega)
    Qgen_new(i) = interp1(OmegaGen,QGen,Omega(i));
end

pos = 0;
neg = 0;

difpos = zeros(length(Omega),1);
difneg = zeros(length(Omega),1);

%Compute the torque difference at each omega
for i=1:length(Omega)
    diff = Q(i)-Qgen_new(i);
    if diff>0
        pos = pos+1;
        difpos(i) = diff;
        difneg(i) = 1e10;
    end
    if diff<0
        neg = neg+1;
        difneg(i) = diff;
        difpos(i) = 1e10;
    end

        
end

%If the torque is higher/lower at each omega, the curves do not interesect
if(pos==0 || neg ==0)
    disp('No convergence at wind speed')
    disp(U)
    Qeq=0;
    Omegaeq=0;
end

if(pos~=0 && neg~=0)
    [M,I] = min(difpos(:))
    [M,Y] = min(abs(difneg(:)))
    i = min(I,Y);
    Delta_Omega = Omega(i+1)-Omega(i);
    m1 = (Q(i+1)-Q(i))/Delta_Omega;
    m2 = (Qgen_new(i+1)-Qgen_new(i))/Delta_Omega;
    Omegaeq = (Qgen_new(i)-Q(i))/(m1-m2);
    Qeq = Omegaeq*m1+Q(i);
end



% if(pos~=0 & neg~=0)
%     for i=1:(length(Omega)-1)
%         Delta_Omega = Omega(i+1)-Omega(i);
%         m1 = (Q(i+1)-Q(i))/Delta_Omega;
%         m2 = (Qgen_new(i+1)-Qgen_new(i))/Delta_Omega;
%         Omegaeq_test = (Qgen_new(i)-Q(i))/(m1-m2);
%         Qeq_test = Omegaeq_test*m1+Q(i);
%         if(Omegaeq_test>Omega(i) & Omegaeq_test<Omega(i+1))
%             Omegaeq = Omegaeq_test;
%             Qeq = Qeq_test;
%         end
%     end
% end

% %Otherwise, find the crossing point
% if(pos~=0 & neg~=0)
%     for i=1:(length(Omega)-1)
%         if((diff(i)>=0 & diff(i+1)<=0)||diff(i)<=0 & diff(i+1)>=0)
%             %Find crossing point
%             Delta_Omega = Omega(i+1)-Omega(i);
%             m1 = (Q(i+1)-Q(i))/Delta_Omega;
%             m2 = (Qgen_new(i+1)-Qgen_new(i))/Delta_Omega;
%             Omegaeq = (Qgen_new(i)-Q(i))/(m1-m2);
%             Qeq = Omegaeq*m1+Q(i);            
%         end
%     end        
% end

 



end

