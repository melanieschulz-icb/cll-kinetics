function [x_out, y_out, z_out] = solution(par, tspanX, tspanY, tspanZ,model)
    params=zeros(1,11);
    params(getRelevantParams(model))=par;
    
    
    d2 = params(1);
    d1 = params(2);
    m = params(3);
    x0 = params(4);
    y0 = params(5);
    c = params(6);
    d3=params(7);
    m2=params(8);
    z0=params(9);
    c2=params(10);
    m3=params(11);

    A = c;
    B = (m+d1+m3);

    A2=c2;
    B2=(m2+d3);

    Cx=A/B;
    CxB1=(x0-Cx);

    Cz=m3*Cx/B2+A2/B2;
    CzB1=m3*CxB1/(B2-B);
    CzB2=z0-Cz-CzB1;

    Cy=(m*Cx+m2*Cz)/d2;
    CyB1=(m*CxB1+m2*CzB1)/(d2-B);
    CyB2=m2*CzB2/(d2-B2);
    Cyd2=y0-Cy-CyB1-CyB2;

    x_out = Cx + CxB1*exp(-B*tspanX);

    z_out=Cz + CzB1*exp(-B*tspanZ)+CzB2*exp(-B2*tspanZ);

    y_out = Cy+CyB1*exp(-B*tspanY)+CyB2*exp(-B2*tspanY)+Cyd2*exp(-d2*tspanY);
% if model == 2
%     [x_out, y_out, z_out]=solutionModel2(par, tspanX, tspanY, tspanZ);
% elseif model == 3
%     [x_out, y_out, z_out]=solutionModel3(par, tspanX, tspanY, tspanZ);
% elseif model == 4
%     [x_out, y_out, z_out]=solutionModel4(par, tspanX, tspanY, tspanZ);
% elseif model == 5
%     [x_out, y_out, z_out]=solutionModel5(par, tspanX, tspanY, tspanZ);
% elseif model == 6
%     [x_out, y_out, z_out]=solutionModel6(par, tspanX, tspanY, tspanZ);
% elseif model == 7
%     [x_out, y_out, z_out]=solutionModel7(par, tspanX, tspanY, tspanZ);
% elseif model == 8
%     [x_out, y_out, z_out]=solutionModel8(par, tspanX, tspanY, tspanZ);
% elseif model == 9
%     [x_out, y_out, z_out]=solutionModel9(par, tspanX, tspanY, tspanZ);
% elseif model == 10
%     [x_out, y_out, z_out]=solutionModel10(par, tspanX, tspanY, tspanZ);
% elseif model == 11
%     [x_out, y_out, z_out]=solutionModel11(par, tspanX, tspanY, tspanZ);
% elseif model == 12
%     [x_out, y_out, z_out]=solutionModel12(par, tspanX, tspanY, tspanZ);
% elseif model == 13
%     [x_out, y_out, z_out]=solutionModel13(par, tspanX, tspanY, tspanZ);
% elseif model == 14
%     [x_out, y_out, z_out]=solutionModel14(par, tspanX, tspanY, tspanZ);
% elseif model == 15
%     [x_out, y_out, z_out]=solutionModel15(par, tspanX, tspanY, tspanZ);
% elseif model == 16
%     [x_out, y_out, z_out]=solutionModel16(par, tspanX, tspanY, tspanZ);
% elseif model == 17
%     [x_out, y_out, z_out]=solutionModel17(par, tspanX, tspanY, tspanZ);
% elseif model == 18
%     [x_out, y_out, z_out]=solutionModel18(par, tspanX, tspanY, tspanZ);
% elseif model == 19
%     [x_out, y_out, z_out]=solutionModel19(par, tspanX, tspanY, tspanZ);
% else
%     print("model not defined");
%     [x_out, y_out, z_out]=[[],[],[]];
end