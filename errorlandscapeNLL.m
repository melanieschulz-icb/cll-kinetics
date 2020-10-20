function stop = errorlandscapeNLL(x,optimValues, state)


stop = false;
switch state
    case 'iter'
        % Find the axes containing the histogram.
     %   NumToNext = ...
     %     findobj(get(gca,'Children'),'Tag','NumberToNextBest');
        
       hold on;
       plot(x(2),optimValues.fval,'*');
    case 'done'
        1
        hold on;
        xlabel('Tissue death rate');
        ylabel('Sum of squared error');
end