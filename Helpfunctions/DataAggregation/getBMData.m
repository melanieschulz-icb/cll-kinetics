function [tspanZAvailable, ZdataAvailable, tspanZ,Zdata,normZdata] = getBMData(patient_nr, startdate, bm_measurements)
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

tspanZAvailable=numPat(1:length(numPat),1)-startdate;
tspanZAvailable(1)=max(tspanZAvailable(1),0);



%extract the measurement values (column 3)
ZdataAvailable=numPat(1:length(numPat),3);

%remove measurements where no data are available (value==nan)
noData=isnan(ZdataAvailable);
tspanZAvailable=tspanZAvailable(~noData);
ZdataAvailable=ZdataAvailable(~noData);
normZdata=ZdataAvailable./ZdataAvailable(1);
ZdataAvailable=ZdataAvailable/100;

tspanZAvailable=tspanZAvailable';
ZdataAvailable=ZdataAvailable';

tspanZ=tspanZAvailable(1:min(bm_measurements,length(tspanZAvailable)));
Zdata=ZdataAvailable(1:length(tspanZ));
normZdata=normZdata';
end