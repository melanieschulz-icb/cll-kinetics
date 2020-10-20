function [dates,Xdata] = getTissueData(patient_nr, cll_volume)
%GETFIRSTPATIENTDATA: gets first patients measurement data after the first
%Ibrutinib dose, reading from the Heavywater-all counts-decoded file
%
%   Detailed explanation goes here
%read the file
[num, txt, ~] = xlsread(strcat(pwd,'/TissueData/',string(patient_nr)));

%find the relevant row, where the tissue values are summed up
[row, ~]=find(~cellfun(@isempty, regexp(txt, 'Total')));

%find the column, where the tissue volume is displayed
[~, col]=find(~cellfun(@isempty, regexp(txt, 'volume')));

%patient9: collumn reads wrong
if patient_nr==9
    col=col+1;
    row=row+1;
end
%extract the numerical values
Xdata=num(row-1,col);
% if patient_nr==28
%     Xdata(1)=1.2*Xdata(1);
% end

%measurement dates
dates=num(:,1);
dates=dates(~isnan(dates));

if patient_nr==4 || patient_nr==9
    dates=dates(2:end);
elseif patient_nr==6 || patient_nr==12 || patient_nr==16 ...
        || patient_nr==23 || patient_nr==24
    dates=[dates; dates(2)];
end
%divide by volume of a cll_cell to estimate the amount of cells
%Xdata=[749.0279;584.3917;491.6138]
Xdata=(Xdata./cll_volume).*10^12;
if patient_nr==6 || patient_nr==7 || patient_nr==12 || patient_nr==16 ...
        || patient_nr==23 || patient_nr==24
    dates=dates(1:2);
    Xdata=Xdata(1:2);
end
end