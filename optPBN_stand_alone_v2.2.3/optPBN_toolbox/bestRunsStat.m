function [cijMeans,cijStds,bestxMean,cijsAll] = bestRunsStat(estim,thresh,threshRelative,SelectedNrRuns)
% [cijMeans,cijStds,bestxMean] = bestRunsStat(estim,thresh,threshRelative,SelectedNrRuns)
% statistics for best runs (bestOptFuncValue threshold)
% select parameter set based on threshold or selected a defined number of parameter set
% calculate mean and standard deviation (SD) of selected parameter sets
%
% Thomas Sauter, University of Luxembourg, 02/2011, thomas.sauter@uni.lu
% Panuwat Trairatphisan, University of Luxembourg 05/2012, panuwat.trairatphisan@uni.lu
% (c) 2014 University of Luxembourg Faculty of Science, Technology and Communication FSTC
% All rights reserved
% GPL version 3.0 to be found at: http://www.gnu.org/licenses/gpl.html

clear history cijsAll
load(estim.fileName,'history','history_old')
history=[history_old; history];

paramNr=estim.paramNr;

if SelectedNrRuns > 0
    Sorted_Result=sortrows(history,-1);
    historySelected=Sorted_Result((end-(SelectedNrRuns-1):end),:);
%     runsThresh=find(historySelected(:,1));
    runsThresh=1:SelectedNrRuns';    
    MeanFVAL=mean(historySelected(:,1))
    NrRunsIncluded=size(historySelected,1)
else
    if threshRelative
        runsThresh=find(history(:,1)<=bestOptFuncValue*thresh)
    else
        runsThresh=find(history(:,1)<=thresh)
    end
    historyThresh=history(runsThresh,:);
    NrRunsIncluded=size(historyThresh,1)
end

if NrRunsIncluded==0
    disp('No runs meeting the defined threshold!')
    disp('Stopping evaluation.')
    cijMeans=[];
    cijStds=[];
    bestxMean=[];
else
    for counter=1:NrRunsIncluded
        if SelectedNrRuns > 0
            cijs=evalCijAndSimulate(estim,[],0,historySelected(counter,(end-estim.paramNr+1):end));
        else
            cijs=evalCijAndSimulate(estim,runsThresh(counter),0);
        end
        cijsAll(:,:,:,counter)=cijs;
    end
    % cijsAll
    for counter=1:size(cijsAll,3)
        disp(['EXPERIMENT' num2str(counter) ':'])
        disp('mean:')
        disp(mean(cijsAll(:,:,counter,:),4))
        disp('std:')
        disp(std(cijsAll(:,:,counter,:),0,4))
    end
    cijMeans=mean(cijsAll,4);
    cijStds=std(cijsAll,0,4);
    
    if SelectedNrRuns > 0
        bestxs=historySelected(runsThresh,(end-estim.paramNr+1):end);
    else
        bestxs=history(runsThresh,(end-estim.paramNr+1):end);
    end
    
    bestxMean=mean(bestxs,1)
end

end

% --- End of script --- %