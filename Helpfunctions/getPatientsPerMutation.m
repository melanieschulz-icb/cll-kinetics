function patients = getPatientsPerMutation(mutationStatus)
%Returns mutated/unmutated/nr patients

if mutationStatus=="mutated"
    patients=[3,6,9,10,14,15,18,21,28:30]; %+30
elseif mutationStatus=="unmutated"
    patients=[1,2,4,5,7,11,13,16,17,19,20,22:26];
elseif mutationStatus=="nr"
    patients=[8,27];
else
    error("specify correct mutation status")
end

