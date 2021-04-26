function stop = errorlandscape(optimValues, state)


stop = false;
switch state
    case 'iter'
        % Find the axes containing the histogram.
     %   NumToNext = ...
     %     findobj(get(gca,'Children'),'Tag','NumberToNextBest');
        
       hold on;
    if length(optimValues.localsolution.X)>0
        plot(optimValues.localsolution.X(2),optimValues.localsolution.Fval,'*');
    end
       
    case 'done'
        hold on;
        xlabel('Tissue death rate');
        ylabel('Sum of squared error');
end