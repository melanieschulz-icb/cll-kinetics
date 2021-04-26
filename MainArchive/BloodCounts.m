[num, txt, ~] = xlsread('Heavywater-all counts-decoded.xlsx','Platelets');
[chars, charDesc,~]=xlsread('Heavy Water patients characterists.xlsx');

[startDate,tspanY,Platelets, HB, Neutroph, RBC] = getBloodCounts(num, txt, chars, charDesc,7);

