function [chars] = getChars(patients)
    allCharsenv=load('Heavy-Water-patients-characterists_IPI_table.mat');
    allChars=allCharsenv.allChars;
    chars=allChars(patients,:);
end

