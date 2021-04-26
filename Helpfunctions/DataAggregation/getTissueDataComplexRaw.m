function [dates,Xdata] = getTissueDataComplexRaw(patient_nr)
%GETFIRSTPATIENTDATA: gets first patients measurement data after the first
%Ibrutinib dose, reading from the Heavywater-all counts-decoded file
%
%   Detailed explanation goes here
%read the file
[num, txt, ~] = xlsread(strcat(pwd,'/TissueData/',string(patient_nr)));

%patient_nr: has a number instead of date -> one row to much
if patient_nr==9
    num=num(2:end,:);
end

%find the relevant row, where the tissue values are summed up
[row, ~]=find(~cellfun(@isempty, regexp(txt, 'Total')));

%find the column, where the tissue volume is displayed
[~, col]=find(~cellfun(@isempty, regexp(txt, 'volume')));

%LNMeasurements

firstLN=num(1:6,8);
secondLN=num(10:15,8);
thirdLN=num(19:24,8);

tooSmall=(firstLN>0);

firstLN=firstLN(tooSmall);
n=length(firstLN);
secondLN=secondLN(tooSmall);
thirdLN=thirdLN(tooSmall);
secondLN(secondLN==0)=1e-3;
thirdLN(thirdLN==0)=1e-3;
%measurement dates
dates=num(:,1);
dates=dates(~isnan(dates));

%patient4: first tissue measurements were done on two different days ->
%chose the second one as only first LN was done different
if patient_nr==4
    dates=dates(2:end);
%patient6,12,16, 23: only two measurements available    
elseif patient_nr==6 || patient_nr==12 || patient_nr==16 ...
        || patient_nr==23 || patient_nr==24
    dates=[dates; dates(2)];
end

%decline rate within the LN
firstRate=(1/(n*(dates(2)-dates(1))))*sum(log(firstLN./secondLN));
secondRate=(1/(n*(dates(3)-dates(2))))*sum(log(secondLN./thirdLN));

%patient9: collumn reads wrong
if patient_nr==9
    col=col+1;
end
%extract the total tissue volumes
Xdata=num(row-1,col);

%compute the declined total volumes, based on the decline rate within the
%LN
Xdata(2)=Xdata(1)*exp(-firstRate*(dates(2)-dates(1)));
Xdata(3)=Xdata(2)*exp(-secondRate*(dates(3)-dates(2)));


%patient6,7,12, 16: only two measurements available 
if patient_nr==6 || patient_nr==7 || patient_nr==12 || patient_nr==16 ...
        || patient_nr==23 || patient_nr==24
    dates=dates(1:2);
    Xdata=Xdata(1:2);
end
end