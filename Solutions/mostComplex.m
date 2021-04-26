function [x_out, y_out, z_out] = mostComplex(params, tspanX, tspanY, tspanZ, model)
% Integrated function of the 3 compartment model with all paths established
%Outputs:
%   time            time span used in the integral
%   x/y/z _out        output from integral for respectively x, y and z
%model model number -> if high number, comartments might be seperated

%Input:
%   par             parameters for integral
%   tspan X/Y/Z       time span for respectively x, y and z
if model <1000
    d2 = params(1);
    d1 = params(2);
    m = params(3);
    x0 = params(4);
    y0 = params(5);
    c = params(6);
    d3=params(7);
    m2=params(8);
    z0=params(9);
    c3=params(10);
    m3=params(11);
    c2=params(12);

    if model <= 154
        A = c;
        B = (m+d1+m3);

        A2=c3;
        B2=(m2+d3);

        Cx=A/B;
        CxB1=(x0-Cx);

        %% difference if one compartment is separated
        if model >= 131 && model <= 154 && d3==0
            Cz=0;
        else
            Cz=m3*Cx/B2+A2/B2;
        end

        %%change::
        if abs(B2-B)<1e-10
            CzB1=m3*CxB1/((B2==B+sign(B2-B))*1e-9);
        else
            CzB1=m3*CxB1/(B2-B);
        end
%       change::: CzB1=m3*CxB1/(B2-B);

        CzB2=z0-Cz-CzB1;
    %%%added c2!!!!!!!!!
    
        Cy=(m*Cx+m2*Cz+c2)/d2;
        %%change:
        if abs(d2-B)<1e-10
            CyB1=(m*CxB1+m2*CzB1)/((d2==B+sign(d2-B))*(1e-9));
        else
            CyB1=(m*CxB1+m2*CzB1)/(d2-B);
        end
        %%change: CyB1=(m*CxB1+m2*CzB1)/(d2-B);
        
        %%change:
        if abs(d2-B2)<1e-10
          CyB2=m2*CzB2/((d2==B2+sign(d2-B2))*(1e-9));
        else
          CyB2=m2*CzB2/(d2-B2);
        end
        %change: CyB2=m2*CzB2/(d2-B2);
        Cyd2=y0-Cy-CyB1-CyB2;

        x_out = Cx + CxB1*exp(-B*tspanX);

        z_out=Cz + CzB1*exp(-B*tspanZ)+CzB2*exp(-B2*tspanZ);

        y_out = Cy+CyB1*exp(-B*tspanY)+CyB2*exp(-B2*tspanY)+Cyd2*exp(-d2*tspanY);

    %% all compartments seperate
    elseif model<=181
        if d1==0
            Cx=0;
        else 
            Cx=c/d1;
        end   
        Cxd1=(x0-Cx);

        if d2==0
            Cy=0;
        else
            Cy=c2/d2;
        end
        Cyd2=(y0-Cy);

        if d3==0
            Cz=0;
        else
            Cz=c3/d3;
        end
        Czd3=(z0-Cz);

        x_out = Cx + Cxd1*exp(-d1*tspanX);

        z_out=Cz + Czd3*exp(-d3*tspanZ);

        y_out = Cy+Cyd2*exp(-d2*tspanY);
    end
else

% Integrated function of the 3 compartment model with all paths established
%Outputs:
%   time            time span used in the integral
%   x/y/z _out        output from integral for respectively x, y and z
%model model number -> if high number, comartments might be seperated

%Input:
%   par             parameters for integral
%   tspan X/Y/Z       time span for respectively x, y and z

    d2 = params(1);
    d1 = params(2);
    m = params(3);
    x0 = params(4);
    y0 = params(5);
    xh = params(6);
    d3=params(7);
    m2=params(8);
    z0=params(9);
    zh=params(10);
    m3=params(11);
    yh=params(12);

    if model <= 2154

        B = (m+d1+m3);
        B2=(m2+d3);

        CxB1=(x0-xh);
        
        if abs(B2-B)<1e-10
            CzB1=m3*CxB1/((B2==B+sign(B2-B))*1e-9);
        else
            CzB1=m3*CxB1/(B2-B);
        end
        CzB2=z0-zh-CzB1; 
        
%         if abs(B2-B)<1e-6
%             z_out= zh + (CzB1*tspanZ+z0-zh).*exp(-B2*tspanZ);   
%         else
%             z_out= zh + CzB1*exp(-B*tspanZ)+CzB2*exp(-B2*tspanZ);            
%         end
        
        if abs(d2-B)<1e-10
            CyB1=(m*CxB1+m2*CzB1)/((d2==B+sign(d2-B))*(1e-9));
        else
            CyB1=(m*CxB1+m2*CzB1)/(d2-B);
        end
        
        if abs(d2-B2)<1e-10
          CyB2=m2*CzB2/((d2==B2+sign(d2-B2))*(1e-9));
        else
          CyB2=m2*CzB2/(d2-B2);
        end
        
%         CzB1=m3*CxB1/(B2-B);
%         CzB2=z0-zh-CzB1;
        
        
%         CyB1=(m*CxB1+m2*CzB1)/(d2-B);

%         
%         CyB2=m2*CzB2/(d2-B2);
        Cyd2=y0-yh-CyB1-CyB2;

        x_out = xh + CxB1*exp(-B*tspanX);

         z_out= zh + CzB1*exp(-B*tspanZ)+CzB2*exp(-B2*tspanZ);

        y_out = yh+CyB1*exp(-B*tspanY)+CyB2*exp(-B2*tspanY)+Cyd2*exp(-d2*tspanY);

    end
end
end