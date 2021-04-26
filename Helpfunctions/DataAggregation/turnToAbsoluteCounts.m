function x_abs = turnToAbsoluteCounts(x_rel,patient)
%turnToAbsoluteCounts: turn relative lymphocyte counts per nanoliter to
%total number of lymphocytes in the body (blood)
%IN:
%   - x_rel: total lymphocyte per nanoliter
%   - patient: patient number
%OUT:
%   - total lymphocyte counts in the blood

volume=estimateBloodVolume(patient);
x_abs=x_rel.*volume;
end