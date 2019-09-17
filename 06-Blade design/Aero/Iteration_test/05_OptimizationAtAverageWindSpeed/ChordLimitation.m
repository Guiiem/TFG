function Chord_New = ChordLimitation(Chord_Old)

len = length(Chord_Old);
Chord_New = zeros(len,1);
Reference_Chord = 0.22;


for i = 1:len
   if Chord_Old(i) <= Reference_Chord
       Chord_New(i) = Chord_Old(i);
   else
       Chord_New(i) = Reference_Chord;
   end    
end



end
