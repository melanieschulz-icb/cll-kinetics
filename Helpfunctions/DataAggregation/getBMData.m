function [tspanZ,Zdata,normZdata] = getBMData(patient_nr, startdate)
%GETFIRSTPATIENTDATA: gets first patients measurement data after the first
%Ibrutinib dose, reading from the Heavywater-all counts-decoded file
%
%   Detailed explanation goes here
%read the file
[num, txt, ~] = xlsread(strcat(pwd,'/Heavywater-BM-CLLtable-no names'));

pat='ACC'+string(patient_nr)+'$';
%extract the relevant data belonging to the patient in question
[idx, idy]=find(~cellfun(@isempty, regexp(txt, pat)));
[accs,~]=find(~cellfun(@isempty, regexp(txt(:,idy:idy+4), 'ACC')));
accsGrPat=accs(accs>idx);
if isempty(accsGrPat)
    endRow=length(num);
else
    endRow=accsGrPat-1;
end
numPat=num(idx:endRow, idy:idy+3);
txtPat=txt(idx:min(endRow, length(txt)), idy:idy+4);

%extract the dates where a measurement was taken (column 1)
%substract the startdate as we want the start date of the therapy to be 0

tspanZ=numPat(1:length(numPat),1)-startdate;
tspanZ(1)=max(tspanZ(1),0);



%extract the measurement values (column 3)
Zdata=numPat(1:length(numPat),3);

%remove measurements where no data are available (value==nan)
noData=isnan(Zdata);
tspanZ=tspanZ(~noData);
Zdata=Zdata(~noData);
normZdata=Zdata./Zdata(1);
Zdata=Zdata/100;

tspanZ=tspanZ';
Zdata=Zdata';
normZdata=normZdata';
end