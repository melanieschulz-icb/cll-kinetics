Xmulti=[];
mutated=[3,6,9,10,14,15,18,21,28:30]; %+30
%mutated is mapped to 1
unmutated=[1,2,4,5,7,11,13,16,17,19,20,22:26];
%unmutated is mapped to 0
nr=[8,27];
%nr is mapped to -1
for patient = [1:11,13:30]
    if patient==30
        [errors, ~,~]=xlsread(strcat(pwd,'/Figures/08_07/ACC29/Errors_ss'));
        [params, ~,~] = xlsread(strcat(pwd,'/Figures/08_07/ACC29/Params_lss'));
    else

        [errors, ~,~]=xlsread(strcat(pwd,'/Figures/08_07/ACC',string(patient),'_01/Errors_ss'));
        [params, ~,~] = xlsread(strcat(pwd,'/Figures/08_07/ACC',string(patient),'_01/Params_lss'));
    end
    bestParams=params(errors(:,2)==min(errors(:,2)),2:end);
    
    if sum(ismember(mutated,patient))>0
        mutStatus=1;
    elseif sum(ismember(unmutated,patient))>0
        mutStatus=0;
    else
        mutStatus=-1;
    end
    Xmulti=[Xmulti; [patient, mutStatus, bestParams]];
        
end

XmultiMut=Xmulti(Xmulti(:,2)==1,:);
XmultiUnmut=Xmulti(Xmulti(:,2)==0,:);
XmultiNr=Xmulti(Xmulti(:,2)==-1,:);

d2=Xmulti(:,3);
d1=Xmulti(:,4);
m=Xmulti(:,5);

x0=Xmulti(:,6);
y0=Xmulti(:,7);
c=Xmulti(:,8);

alph=d1+m;
Cx=(d1./(m+d1)).*c;
z=(m./alph).*(1+(Cx./x0).*(alph.*(1./(m+d1)+1./d2)-1));

writetable(array2table([Xmulti,alph,z], 'VariableNames',{'ACC','status','d2', 'd1','m','x0','y0','c','alph','z'}),'bestParams_lss.xlsx');