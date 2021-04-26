function stop = errorlandscapeLSQ(x,optimValues, state)

hold on;
stop = false;
switch state
    case 'iter'
        % Find the axes containing the histogram.
     %   NumToNext = ...
     %     findobj(get(gca,'Children'),'Tag','NumberToNextBest');
        
       hold on;
       plot(x(2),optimValues.resnorm,'*');
       hold on;
%    case 'done'
%        1
%        hold on;
%        xlabel('Tissue death rate');
%        ylabel('Sum of squared error');
end