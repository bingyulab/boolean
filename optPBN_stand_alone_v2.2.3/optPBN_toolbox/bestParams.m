function [bestOptFuncValue,bestRun,bestx] = bestParams(estim,plotFLAG)
% [bestOptFuncValue,bestRun,bestx] = bestParams(estim,plotFLAG)
% identification history: identify bast parameter set & plot (optional)
%
% Thomas Sauter, University of Luxembourg, 02/2011, thomas.sauter@uni.lu
% Panuwat Trairatphisan, University of Luxembourg, 04/2014, panuwat.trairatphisan@uni.lu
% (c) 2014 University of Luxembourg Faculty of Science, Technology and Communication FSTC
% All rights reserved
% GPL version 3.0 to be found at: http://www.gnu.org/licenses/gpl.html

load(estim.fileName,'history','history_old')
history=[history_old; history];

[bestOptFuncValue,bestRun]=(min(history(:,1)));
disp('best optimal function value, pdf_ic_mean & parameter:')
disp(history(bestRun,:))

if plotFLAG==1
    figure, plot(history(:,1)), title('FVAL during optimisation')
    figure, plot(sort(history(:,1),'descend')), title('FVAL all sorted')
end

bestx=history(bestRun,(end-estim.paramNr+1):end);

end

% --- End of script --- %