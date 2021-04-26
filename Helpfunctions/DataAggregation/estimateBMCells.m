function [bm_cells] = estimateBMCells(patient)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
chars=getChars(patient);
if chars.SEX=="M"
    factor=7.8;
else
    factor=6.7;
end
bm_cells=factor*chars.Weight_kg_*1e9;
end

