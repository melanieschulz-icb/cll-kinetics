function [startDate,tspanY,Ydata] = getPatientData(num, txt,patient_nr, beforeTherapy)
%GETFIRSTPATIENTDATA: gets first patients measurement data after the first
%Ibrutinib dose or before therapy, reading from the Heavywater-all counts-decoded file
%

if ~exist('beforeTherapy','var')
    beforeTherapy=false;
end

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

% find the index of the therapy or HW start, by using regexp

[idxHW, ~]=find(~cellfun(@isempty, regexp(txtPat,'HW ?'))); 
[idxIb, ~]=find(~cellfun(@isempty, regexp(txtPat,'brutinib ?')));

%needs to be shiftet as headers are missing in the numbers table
idxHW=idxHW-2;
idxIb=idxIb-2;

%idx is the row where the interesting values begin
if beforeTherapy
    idx=idxHW;
else
    idx=idxIb;
end

%extract the dates where a measurement was taken (column 4)
%substract the first value as we want the start date of the therapy to be 0
tspanY=numPat(idx:length(numPat),4)-numPat(idxIb,4);
%extract the measurement values (column 2)
Ydata=numPat(idx:length(numPat),2);
%startDate (column 3)
startDate=numPat(idxIb,3);
%remove measurements where no data are available (value==nan)
noData=isnan(Ydata);
tspanY=tspanY(~noData);


% Ydata=Ydata(~noData);

Ydata=turnToAbsoluteCounts(Ydata(~noData),patient_nr);
end

