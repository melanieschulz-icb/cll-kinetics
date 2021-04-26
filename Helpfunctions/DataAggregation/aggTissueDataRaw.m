function [tspanXAvailable, XdataAvailable, tspanX, xdata]= aggTissueDataRaw(patient, use_complex_tss, startDate, tss_measurements)        
            if use_complex_tss==1    
                [tspanXAvailable,XdataAvailable]=getTissueDataComplexRaw(patient);
            else
                [tspanXAvailable, XdataAvailable]=getTissueDataRaw(patient);
            end
            
            tspanXAvailable=(tspanXAvailable-startDate)';
            tspanXAvailable(tspanXAvailable<0)=0;
            tspanX=tspanXAvailable(1:min(tss_measurements,length(tspanXAvailable)));
            xdata=XdataAvailable(1:length(tspanX))';
            
            if patient==24 && tss_measurements==2
                %second measurement missing -> if we ant to have measurements ...
                %on the first and second timepoint, we have to remove the second datapoint...
                %as it belongs to the third measurement!!!
                tspanX=tspanX(1:1);
                xdata=xdata(1:1);
            end
end