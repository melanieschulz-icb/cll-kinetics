function [startDate, tspanY, ydata, ty_outlier, y_outlier] = aggBloodData(num, txt,patient, beforeTherapy)
    if ~exist('beforeTherapy','var')
        beforeTherapy=false;
    end
        
            [startDate,tspanY, ydata]=getPatientData(num, txt,patient, beforeTherapy);
            tspanY=tspanY';
            ydata=ydata';
            ty_outlier=[];
            y_outlier=[];
            
            %fix data
            if patient == 3
                outlier=4;
                ty_outlier=tspanY(outlier);
                tspanY=[tspanY(1:outlier-1),tspanY(outlier+1: end)];
                y_outlier=ydata(outlier);
                ydata=[ydata(1:outlier-1),ydata(outlier+1: end)];
            elseif patient == 7
                %datum falsch eingetragen (+1Jahr)
                tspanY(5)=tspanY(5)-365;
            elseif patient == 8
                outlier=4;
                ty_outlier=tspanY(outlier);
                tspanY=[tspanY(1:outlier-1),tspanY(outlier+1: end)];
                y_outlier=ydata(outlier);
                ydata=[ydata(1:outlier-1),ydata(outlier+1: end)];
            end
end
            