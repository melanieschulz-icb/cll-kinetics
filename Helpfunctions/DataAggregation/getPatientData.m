function [startDate,tspanY,Ydata] = getPatientData(num, txt, chars, charDescs,patient_nr)
%GETFIRSTPATIENTDATA: gets first patients measurement data after the first
%Ibrutinib dose, reading from the Heavywater-all counts-decoded file
%
%   Detailed explanation goes here
pat='ACC '+string(patient_nr)+'$';
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

% find the index of the therapy start, by using regexp
[idx, ~]=find(~cellfun(@isempty, regexp(txtPat,'brutinib ?')));
%idx is the row where the interesting values begin
%needs to be shiftet as headers are missing in the numbers table
idx=idx-2;
%extract the dates where a measurement was taken (column 4)
%substract the first value as we want the start date of the therapy to be 0
tspanY=numPat(idx:length(numPat),4)-numPat(idx,4);
%extract the measurement values (column 2)
Ydata=numPat(idx:length(numPat),2);
%startDate (column 3)
startDate=numPat(idx,3);
%remove measurements where no data are available (value==nan)
noData=isnan(Ydata);
tspanY=tspanY(~noData);

%get factors for bloodVolume
weight=chars(patient_nr,7);
sex=charDescs(patient_nr,5);

if string(sex)=='M'
    factor=70;
else
    factor=61;
end


Ydata=factor*weight*10^6*Ydata(~noData);

end

