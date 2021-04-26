function volume = estimateBloodVolume(patient)
%turnToAbsoluteCounts: turn relative lymphocyte counts per nanoliter to
%total number of lymphocytes in the body (blood)
%IN:
%   - x_rel: total lymphocyte per nanoliter
%   - patient: patient number
%OUT:
%   - total lymphocyte counts in the blood

chars=getChars(patient);

sex=chars.SEX;
factor=61*ones(length(sex),1);

for i = 1:length(sex)
    if sex{i}=='M'
        factor(i)=70;
    end
end

weight=chars.Weight_kg_;

volume=factor.*weight.*1e6;
end