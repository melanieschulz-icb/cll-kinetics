function bm_totalCells = getBMTotalCells(cll_cells,patients)
%getBMTotalCells Get the number of total cells within the BM at timepoint
%0, having the number of cll cells given
%   IN: 
%       - cell_cells: absolute number of tumor cells at timepoint 0
%       - patients: patients, to determine the percentage of bm cells that
%           are tumor cells
percentage=ones(length(patients),1);
for i=1:length(patients)
    [~, ~, ~,Zdata,~] = getBMData(patients(i), 0, 1);
    percentage(i)=Zdata;
end
bm_totalCells=cll_cells./Zdata;
end

