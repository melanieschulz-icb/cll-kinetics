function [tspanXAvailable, XdataAvailable, tspanX, xdata]= aggTissueData(patient, cll_volume, use_complex_tss, startDate, tss_measurements)        
            if use_complex_tss==1    
                [tspanXAvailable,XdataAvailable]=getTissueDataComplex(patient,cll_volume);
            else
                [tspanXAvailable, XdataAvailable]=getTissueData(patient, cll_volume);
            end
            
            tspanXAvailable=(tspanXAvailable-startDate)';
            tspanXAvailable(tspanXAvailable<0)=0;
            tspanX=tspanXAvailable(1:tss_measurements);
            xdata=XdataAvailable(1:tss_measurements)';
            
            if patient==24
                %second measurement missing!!!
                tspanX=tspanX(1:1);
                xdata=xdata(1:1);
            end
end