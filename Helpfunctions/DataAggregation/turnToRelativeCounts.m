function x_rel = turnToRelativeCounts(x_abs,patient)
%turnToRelativeCounts: turn total lymphocyte counts to counts per
%nanoliter
%IN:
%   - x_abs: total lymphocyte counts in the body
%   - patient: patient number
%OUT:
%   - lymphocyte counts per nanoliter

volume=estimateBloodVolume(patient);

x_rel=x_abs./volume;
end

