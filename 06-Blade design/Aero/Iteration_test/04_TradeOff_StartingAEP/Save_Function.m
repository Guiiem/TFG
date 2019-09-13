function [] = Save_Function(Operation,Blade,Ind,Data)

% Save_Flag = 1; %Select whether we want to save 
% Save_Name = 'Results_1';
% 
% Save_Name_CT = strcat(Save_Name,'Chord_Twist','.txt');
% 
% fileID = fopen(Save_Name_CT,'w');
% fprintf(fileID,'%6s %6s \n','Chord','Twist');
% for i=1:Data.no_genes
%     fprintf(fileID,'%5.4f %5.4f \n', Blade.Chord(Ind,i),Blade.Twist(Ind,i));
% end
% fclose(fileID);

save('C:\Users\guill\Dropbox\TFG\06-Blade design\Aero\Iteration_test\!Results\1_EOGEN_500W_1.2m_Itrt1.mat')


end