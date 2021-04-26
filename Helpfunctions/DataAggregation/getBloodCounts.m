function [startDate,tspanY,Platelets,PlateletsScaled, HB,HBScaled, Neutroph,NeutrophScaled, RBC, RBCScaled] = getBloodCounts(num, txt, chars, charDescs,patient_nr, beforeTherapy)
%GETFIRSTPATIENTDATA: gets first patients measurement data after the first
%Ibrutinib or HW dose, reading from the Heavywater-all counts-decoded file
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
numPat=num(idx:endRow, idy:idy+5);
txtPat=txt(idx:min(endRow, length(txt)), idy:idy+6);

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
tspanY=numPat(idx:length(numPat),5)-numPat(idxIb,5);
%extract the measurement values (column 2)
Platelets=numPat(idx:length(numPat),1);
HB=numPat(idx:length(numPat),2);
Neutroph=numPat(idx:length(numPat),3);
RBC=numPat(idx:length(numPat),4);
%startDate (column 3)
startDate=numPat(idxIb,5);
%remove measurements where no data are available (value==nan)
noData=isnan(Platelets);
tspanY=tspanY(~noData);

volume=estimateBloodVolume(patient_nr);

Platelets=Platelets(~noData);
HB=HB(~noData);
Neutroph=Neutroph(~noData);
RBC=RBC(~noData);

PlateletsScaled=volume*Platelets;
HBScaled=volume*HB;
NeutrophScaled=volume*Neutroph;
RBCScaled=volume*RBC;



end