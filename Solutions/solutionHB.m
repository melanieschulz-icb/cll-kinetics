function [hb_out] = solutionHB(par, tspanYPlatelets ,bloodmodel, model)
    
    if isempty(model)
        model=1;
    end
    
    %put params into modelParams
    %HB_0, HB_max, k, d
    if model == 1
        HB_0=par(1);
        HB_max=par(2);
        k=par(3);

        hb_out=HB_max-(HB_max-HB_0)*exp(-k*tspanYPlatelets);
    elseif model==2
        HB_0=par(1);
        HB_max=par(2);
        k=par(3);
        d=par(4);
        hb_out=HB_max-(HB_max-HB_0)*exp(d*tspanYPlatelets-k/2*tspanYPlatelets.^2);
    elseif model==3
        bloodParams=zeros(1,12);
        bloodParams(getRelevantParams(bloodmodel))=par(1:end-4);
        

        d1 = bloodParams(2);
        m = bloodParams(3);
        x0 = bloodParams(4);

        xh = bloodParams(6);
        d3=bloodParams(7);
        m2=bloodParams(8);
        z0=bloodParams(9);
        zh=bloodParams(10);
        m3=bloodParams(11);

        
        HB_0=par(end-3);
        HB_max=par(end-2);
        k=par(end-1);
        d=par(end);
        
        B = (m+d1+m3);
        B2=(m2+d3);

        CxB1=(x0-xh);

        if B2-B==0
            B2=B2+0.0001;
        end
        
        CzB1=m3*CxB1/(B2-B);
        CzB2=z0-zh-CzB1;
        
        z_out= zh + CzB1*exp(-B*tspanYPlatelets)+CzB2*exp(-B2*tspanYPlatelets);
        z_scaled=z_out/z0;

        hb_out=HB_max-(HB_max-HB_0)*exp(k*z_scaled-(k-d)*tspanYPlatelets-k);
    elseif model==4
        bloodParams=zeros(1,12);
        bloodParams(getRelevantParams(bloodmodel))=par(1:end-3);
        

        d1 = bloodParams(2);
        m = bloodParams(3);
        x0 = bloodParams(4);

        xh = bloodParams(6);
        d3=bloodParams(7);
        m2=bloodParams(8);
        z0=bloodParams(9);
        zh=bloodParams(10);
        m3=bloodParams(11);

        
        HB_0=par(end-2);
        HB_max=par(end-1);
        k=par(end);

        
        B = (m+d1+m3);
        B2=(m2+d3);

        CxB1=(x0-xh);

        if B2-B==0
            B2=B2+0.0001;
        end
        
        CzB1=m3*CxB1/(B2-B);
        CzB2=z0-zh-CzB1;
        
        z_out= zh + CzB1*exp(-B*tspanYPlatelets)+CzB2*exp(-B2*tspanYPlatelets);
        z_scaled=z_out/z0;

        hb_out=HB_max-(HB_max-HB_0)*exp(k*z_scaled-k*tspanYPlatelets-k);
  
    elseif model==5
        bloodParams=zeros(1,12);
        bloodParams(getRelevantParams(bloodmodel))=par(1:end-3);
        

        d1 = bloodParams(2);
        m = bloodParams(3);
        x0 = bloodParams(4);

        xh = bloodParams(6);
        d3=bloodParams(7);
        m2=bloodParams(8);
        z0=bloodParams(9);
        zh=bloodParams(10);
        m3=bloodParams(11);

        
        HB_0=par(end-2);
        HB_max=par(end-1);
        k=par(end);

        
%         B = (m+d1+m3);
%         B2=(m2+d3);
% 
%         CxB1=(x0-xh);
% 
%         if B2-B==0
%             B2=B2+0.0001;
%         end
%         
%         CzB1=m3*CxB1/(B2-B);
%         CzB2=z0-zh-CzB1;
        
%         z_out= zh + CzB1*exp(-B*tspanYPlatelets)+CzB2*exp(-B2*tspanYPlatelets);
        [~,~,z_out_test]=solution(par(1:end-3),[], [],tspanYPlatelets, bloodmodel);
        z_scaled=z_out_test/z0;
%         z_scaled=z_out/z0;

        hb_out=min(HB_max-(HB_max-HB_0)*exp(k*z_scaled-k),100000);
        
        
     elseif model==6
        bloodParams=zeros(1,12);
        bloodParams(getRelevantParams(bloodmodel))=par(1:end-4);
        

        d1 = bloodParams(2);
        m = bloodParams(3);
        x0 = bloodParams(4);

        xh = bloodParams(6);
        d3=bloodParams(7);
        m2=bloodParams(8);
        z0=bloodParams(9);
        zh=bloodParams(10);
        m3=bloodParams(11);

        
        HB_0=par(end-3);
        HB_max=par(end-2);
        k=par(end-1);
        d=par(end);
        
        B = (m+d1+m3);
        B2=(m2+d3);

        CxB1=(x0-xh);

        if B2-B==0
            B2=B2+0.0001;
        end
        
        CzB1=m3*CxB1/(B2-B);
        CzB2=z0-zh-CzB1;
        
        z_out= zh + CzB1*exp(-B*tspanYPlatelets)+CzB2*exp(-B2*tspanYPlatelets);
        z_scaled=z_out/z0;
       % z_scaled=z_out;

        hb_out=HB_max-(HB_max-HB_0)*exp(k*z_scaled+d*tspanYPlatelets-k);
        
     elseif model==7
        bloodParams=zeros(1,12);
        bloodParams(getRelevantParams(bloodmodel))=par(1:end-3);
        

        d1 = bloodParams(2);
        m = bloodParams(3);
        x0 = bloodParams(4);

        xh = bloodParams(6);
        d3=bloodParams(7);
        m2=bloodParams(8);
        z0=bloodParams(9);
        zh=bloodParams(10);
        m3=bloodParams(11);

        
        HB_0=par(end-2);
    
        k=par(end-1);
        d=par(end);
        
        B = (m+d1+m3);
        B2=(m2+d3);

        CxB1=(x0-xh);

        if B2-B==0
            B2=B2+0.0001;
        end
        
        CzB1=m3*CxB1/(B2-B);
        CzB2=z0-zh-CzB1;
        
        z_out= zh + CzB1*exp(-B*tspanYPlatelets)+CzB2*exp(-B2*tspanYPlatelets);
        z_scaled=z_out/z0;
       % z_scaled=z_out;

        hb_out=HB_0+(k*(1-z_scaled)-d).*tspanYPlatelets;
        
    elseif model==8
        
        bloodParams=zeros(1,12);
        bloodParams(getRelevantParams(bloodmodel))=par(1:end-3);
        

%        d2 = bloodParams(1);
        d1 = bloodParams(2);
        m = bloodParams(3);
        x0 = bloodParams(4);
%        y0 = bloodParams(5);
        xh = bloodParams(6);
        d3=bloodParams(7);
        m2=bloodParams(8);
        z0=bloodParams(9);
        zh=bloodParams(10);
        m3=bloodParams(11);
%        yh=bloodParams(12);

        
        HB_0=par(end-2);
%         HB_max=par(end-2);
        k=par(end-1);
        d=par(end);
        
        B = (m+d1+m3);
        B2=(m2+d3);

        CxB1=(x0-xh);

        if B2-B==0
            B2=B2+0.0001;
        end
        
        CzB1=m3*CxB1/(B2-B);
        CzB2=z0-zh-CzB1;
        
        hb_h=k*(1-zh/z0)/d;
        ChbB1=k*CzB1/(z0*(d-B));
        ChbB2=k*CzB2/(z0*(d-B2));
        Chbd=HB_0-hb_h+ChbB1+ChbB2;
        
        hb_out=hb_h-ChbB1*exp(-B*tspanYPlatelets)-ChbB2*exp(-B2*tspanYPlatelets)+Chbd*exp(-d*tspanYPlatelets);
        

    elseif model==9
        bloodParams=zeros(1,12);
        bloodParams(getRelevantParams(bloodmodel))=par(1:end-4);
        

%        d2 = bloodParams(1);
        d1 = bloodParams(2);
        m = bloodParams(3);
        x0 = bloodParams(4);
%        y0 = bloodParams(5);
        xh = bloodParams(6);
        d3=bloodParams(7);
        m2=bloodParams(8);
        z0=bloodParams(9);
        zh=bloodParams(10);
        m3=bloodParams(11);
%        yh=bloodParams(12);

        
        HB_0=par(end-3);
        HB_max=par(end-2);
        k=par(end-1);
        d=par(end);
        
        B = (m+d1+m3);
        B2=(m2+d3);

        CxB1=(x0-xh);

        if B2-B==0
            B2=B2+0.0001;
        end
        
        CzB1=m3*CxB1/(B2-B);
        CzB2=z0-zh-CzB1;
        
        hb_h=k*(zh/z0-1)/d + HB_max;
        
        ChbB1=k*CzB1/(z0*(d+B));
        ChbB2=k*CzB2/(z0*(d+B2));
        Chbd=HB_0-(hb_h+ChbB1+ChbB2);
        
        hb_out=hb_h+ChbB1*exp(-B*tspanYPlatelets)+ChbB2*exp(-B2*tspanYPlatelets)+Chbd*exp(d*tspanYPlatelets);
   
    elseif model==10
        bloodParams=zeros(1,12);
        bloodParams(getRelevantParams(bloodmodel))=par(1:end-4);
        

%        d2 = bloodParams(1);
        d1 = bloodParams(2);
        m = bloodParams(3);
        x0 = bloodParams(4);
%        y0 = bloodParams(5);
        xh = bloodParams(6);
        d3=bloodParams(7);
        m2=bloodParams(8);
        z0=bloodParams(9);
        zh=bloodParams(10);
        m3=bloodParams(11);
%        yh=bloodParams(12);

        
        HB_0=par(end-3);
        HB_max=par(end-2);
        k=par(end-1);
        d=par(end);
        
        B = (m+d1+m3);
        B2=(m2+d3);

        CxB1=(x0-xh);

        if B2-B==0
            B2=B2+0.0001;
        end
        
        CzB1=m3*CxB1/(B2-B);
        CzB2=z0-zh-CzB1;
        
        fun = @(x) (1./(HB_0.*exp((k/z0).*(CzB1/B+CzB2/B2)+(d-k+(k/z0)*zh).*x-(k/z0).*((CzB1/B).*exp(-B.*x)+(CzB2/B2).*exp(-B2.*x))))).*k.*(1-zh/z0-(CzB1/z0).*exp(-B.*x)-(CzB2/z0).*exp(-B2.*x)).*HB_max;
        
        int=zeros(1,length(tspanYPlatelets));
        counter=1;
        for t=tspanYPlatelets
            int(counter)=integral(fun,0,t);
            counter=counter+1;
        end
        
        hb_out=HB_0.*exp(k/z0.*(CzB1/B+CzB2/B2)+(d-k+k/z0*zh).*tspanYPlatelets-...
            k/z0.*(CzB1/B*exp(-B.*tspanYPlatelets)+CzB2/B2.*exp(-B2*tspanYPlatelets))).*...
            (HB_0+int);
        
        
        
    elseif model==11
        HB_0=par(1);
        HB_max=par(2);
        k=par(3);
        d=par(4);
        hb_out=HB_max-(HB_max-HB_0)*exp(-d*tspanYPlatelets+k/2*tspanYPlatelets.^2);
    end
end