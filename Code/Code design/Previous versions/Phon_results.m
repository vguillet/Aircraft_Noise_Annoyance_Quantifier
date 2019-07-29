% compare Aspl to Phonlines
phonresult = zeros(29,1);
 for i = 5:29
     for j = 1:900
         
     
         if aspl(i) <= phontable(j,i)
             phonresult(i,1) = j/10;
             break
         end
        
     end
 end
 
         