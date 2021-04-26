function [x_out, y_out, z_out] = solution(par, tspanX, tspanY, tspanZ,model)
 
        params=zeros(1,12);
    if model==0
        z_out=[];
        [x_out, y_out]=intLea(par,tspanX,tspanY);
    elseif model==0.5
        z_out=[];
        [x_out, y_out]=intBurger(par,tspanX,tspanY);
    elseif model<=2
        z_out=[];
        params=zeros(1,7);
        params(getRelevantParams(model))=par;
        params=params([1:6,12]);
        [x_out, y_out]=intLeaSeparation(params,tspanX,tspanY);
    elseif model <= 181 
        %LT->PB
        %LT->BM
        %BM->PB
         params(getRelevantParams(model))=par;
        [x_out,y_out, z_out]=mostComplex(params, tspanX, tspanY, tspanZ, model);

    elseif model<=309 
        %BM->PB
        %BM->LT
        %LT->PB
        %-> swap x,z
        params(getRelevantParams(model-179))=par;
        [z_out,y_out, x_out]=mostComplex(params, tspanZ, tspanY, tspanX,model-179);
    elseif model <=437 
        %LT->BM
        %LT->PB
        %PB->BM
        %->swap y,z
        params(getRelevantParams(model-307))=par;
        [x_out,z_out, y_out]=mostComplex(params, tspanX, tspanZ, tspanY, model-307);
    elseif model <=565
        %BM->LT
        %BM->PB
        %PB->LT
        %z->x->y
        params(getRelevantParams(model-435))=par;
        [z_out,x_out, y_out]=mostComplex(params, tspanZ, tspanX, tspanY,model-435);
    elseif model <=693
        %PB->BM
        %PB->LT
        %LT->BM
        %y->x->z
        params(getRelevantParams(model-563))=par;
        [y_out,z_out, x_out]=mostComplex(params, tspanY, tspanZ, tspanX, model-563);
    elseif model <=821
        %PB->LT
        %PB->BM
        %BM->LT
        %swap x,y
        params(getRelevantParams(model-691))=par;
        [y_out,x_out, z_out]=mostComplex(params, tspanY, tspanX, tspanZ,model-691);
        
    elseif model <=1021
        %LT->PB
        %LT->BM
        %BM->PB
         params(getRelevantParams(model))=par;
        [x_out,y_out, z_out]=mostComplex(params, tspanX, tspanY, tspanZ, model);
    elseif model<=1037
        %BM->PB
        %BM->LT
        %LT->PB
        %-> swap x,z
        params(getRelevantParams(model-19))=par;
        [z_out,y_out, x_out]=mostComplex(params, tspanZ, tspanY, tspanX,model-19);
    elseif model <=1053 
        %LT->BM
        %LT->PB
        %PB->BM
        %->swap y,z
        params(getRelevantParams(model-35))=par;
        [x_out,z_out, y_out]=mostComplex(params, tspanX, tspanZ, tspanY, model-35);
    elseif model <=1069
        %BM->LT
        %BM->PB
        %PB->LT
        %z->x->y
        params(getRelevantParams(model-51))=par;
        [z_out,x_out, y_out]=mostComplex(params, tspanZ, tspanX, tspanY,model-51);
    elseif model <=1085
        %PB->BM
        %PB->LT
        %LT->BM
        %y->x->z
        params(getRelevantParams(model-67))=par;
        [y_out,z_out, x_out]=mostComplex(params, tspanY, tspanZ, tspanX, model-67);
    elseif model <=1101
        %PB->LT
        %PB->BM
        %BM->LT
        %swap x,y
        params(getRelevantParams(model-83))=par;
        [y_out,x_out, z_out]=mostComplex(params, tspanY, tspanX, tspanZ,model-83);
    elseif model <=2021
        %LT->PB
        %LT->BM
        %BM->PB
         params(getRelevantParams(model))=par;
        [x_out,y_out, z_out]=mostComplex(params, tspanX, tspanY, tspanZ, model);
    elseif model<=2037
        %BM->PB
        %BM->LT
        %LT->PB
        %-> swap x,z
        %params(getRelevantParams(model-19))=par;
        %[z_out,y_out, x_out]=mostComplex(params, tspanZ, tspanY, tspanX,model-19);
        pos_with_zh=getRelevantParams(model-1019);
        names=getVarNames(model-1000);
        ind_zh=find(names(getRelevantParams(model-1000))=="BM_healthy");
        params(pos_with_zh([1:end]~=ind_zh))=par;
%        params([1     2     3     4     5     7     9  10  11    12])=par;
        [z_out,y_out, x_out]=mostComplex(params, tspanZ, tspanY, tspanX,model-19);
    elseif model <=2053 
        %LT->BM
        %LT->PB
        %PB->BM
        %->swap y,z
        pos_with_zh=getRelevantParams(model-1035);
        names=getVarNames(model-1000);
        ind_zh=find(names(getRelevantParams(model-1000))=="BM_healthy");
        params(pos_with_zh([1:end]~=ind_zh))=par;
        [x_out,z_out, y_out]=mostComplex(params, tspanX, tspanZ, tspanY, model-35);
    elseif model <=2069
        %BM->LT
        %BM->PB
        %PB->LT
        %z->x->y
        pos_with_zh=getRelevantParams(model-1051);
        names=getVarNames(model-1000);
        ind_zh=find(names(getRelevantParams(model-1000))=="BM_healthy");
        params(pos_with_zh([1:end]~=ind_zh))=par;
        [z_out,x_out, y_out]=mostComplex(params, tspanZ, tspanX, tspanY,model-51);
    elseif model <=2085
        %PB->BM
        %PB->LT
        %LT->BM
        %y->x->z
        pos_with_zh=getRelevantParams(model-1067);
        names=getVarNames(model-1000);
        ind_zh=find(names(getRelevantParams(model-1000))=="BM_healthy");
        params(pos_with_zh([1:end]~=ind_zh))=par;
        [y_out,z_out, x_out]=mostComplex(params, tspanY, tspanZ, tspanX, model-67);
    elseif model <=2101
        %PB->LT
        %PB->BM
        %BM->LT
        %swap x,y
        pos_with_zh=getRelevantParams(model-1083);
        names=getVarNames(model-1000);
        ind_zh=find(names(getRelevantParams(model-1000))=="BM_healthy");
        params(pos_with_zh([1:end]~=ind_zh))=par;
        [y_out,x_out, z_out]=mostComplex(params, tspanY, tspanX, tspanZ,model-83);
        
%     else
%         params=zeros(1,12);
%     
% %         if model < 100
% %             params(getRelevantParams(model))=par;
% %         else
% %             params(getRelevantParams(model-100))=par;
% %         end
% 
%         params(getRelevantParams(model))=par;
%         
%         d2 = params(1);
%         d1 = params(2);
%         m = params(3);
%         x0 = params(4);
%         y0 = params(5);
%         c = params(6);
%         d3=params(7);
%         m2=params(8);
%         z0=params(9);
%         c3=params(10);
%         m3=params(11);
%         c2=params(12);
% 
%         if model <= 154
%             A = c;
%             B = (m+d1+m3);
%             
%             A2=c3;
%             B2=(m2+d3);
% 
%             Cx=A/B;
%             CxB1=(x0-Cx);
%             
%             %% difference if one compartment is seperated
%             if model >= 131 && model <= 154 && d3==0
%                 Cz=0;
%             else
%                 Cz=m3*Cx/B2+A2/B2;
%             end
%             
%             %%
%             CzB1=m3*CxB1/(B2-B);
%             CzB2=z0-Cz-CzB1;
% %%%added c2!!!!!!!!!
%             Cy=(m*Cx+m2*Cz+c2)/d2;
%             CyB1=(m*CxB1+m2*CzB1)/(d2-B);
%             CyB2=m2*CzB2/(d2-B2);
%             Cyd2=y0-Cy-CyB1-CyB2;
% 
%             x_out = Cx + CxB1*exp(-B*tspanX);
% 
%             z_out=Cz + CzB1*exp(-B*tspanZ)+CzB2*exp(-B2*tspanZ);
% 
%             y_out = Cy+CyB1*exp(-B*tspanY)+CyB2*exp(-B2*tspanY)+Cyd2*exp(-d2*tspanY);
%             
%         %% all compartments seperate
%         elseif model<=181
%             if d1==0
%                 Cx=0;
%             else 
%                 Cx=c/d1;
%             end   
%             Cxd1=(x0-Cx);
%             
%             if d2==0
%                 Cy=0;
%             else
%                 Cy=c2/d2;
%             end
%             Cyd2=(y0-Cy);
%             
%             if d3==0
%                 Cz=0;
%             else
%                 Cz=c3/d3;
%             end
%             Czd3=(z0-Cz);
%             
%             x_out = Cx + Cxd1*exp(-d1*tspanX);
% 
%             z_out=Cz + Czd3*exp(-d3*tspanZ);
% 
%             y_out = Cy+Cyd2*exp(-d2*tspanY);
%         %%    
%         end
%         else
%     %        m3 from bm to lt -> interchange z and x 
%             A = c2;
%             B = (m2+d3+m3);
% 
%             A2=c;
%             B2=(m+d1);
% 
%             Cz=A/B;
%             CzB1=(z0-Cz);
% 
%             Cx=m3*Cz/B2+A2/B2;
%             CxB1=m3*CzB1/(B2-B);
%             CxB2=x0-Cx-CxB1;
% 
%             Cy=(m2*Cz+m*Cx)/d2;
%             CyB1=(m*CxB1+m2*CzB1)/(d2-B);
%             CyB2=m*CxB2/(d2-B2);
%             Cyd2=y0-Cy-CyB1-CyB2;
% 
%             z_out = Cz + CzB1*exp(-B*tspanZ);
% 
%             x_out=Cx + CxB1*exp(-B*tspanX)+CxB2*exp(-B2*tspanX);
% 
%             y_out = Cy+CyB1*exp(-B*tspanY)+CyB2*exp(-B2*tspanY)+Cyd2*exp(-d2*tspanY);
%         end
    end
end