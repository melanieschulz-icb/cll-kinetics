function [x_out, y_out] = intLeaSeparation(par, tspanX, tspanY)
%     params=zeros(1,10);
%     params(getRelevantParams(model,stabilization))=par;
    params=par;    

    d2 = params(1);
    d1 = params(2);
    m = params(3);
    x0 = params(4);
    y0 = params(5);
    xh = params(6);
    yh = params(7);
    
    CxB1=(x0-xh);
    B = m+d1;
    
    if abs(d2-B)<1e-10
        CyB1=m*CxB1/((d2==B+sign(d2-B))*1e-9);
    else
        CyB1=m*CxB1/(d2-B);
    end
    
    Cyd2=y0-yh-CyB1;
    
    x_out = xh + CxB1*exp(-B*tspanX);

    y_out= yh + CyB1*exp(-B*tspanY)+Cyd2*exp(-d2*tspanY);    



%         if abs(d2-B)<1e-10
%             CY=(m*x0)/((d2==B+sign(d2-B))*(1e-9));
%         else
%             CY = (m*x0)/(-B+d2);
%         end
%         
%         
% 
%         x_out = xh + x0*exp(-B*tspanX);
% 
%         y_out = yh + CY*exp(-B*tspanY) + (y0-CY)*exp(-d2*tspanY);
    
end