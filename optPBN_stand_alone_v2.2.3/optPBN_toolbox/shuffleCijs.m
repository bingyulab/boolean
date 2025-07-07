function [cijMeansShuff,cijStdsShuff] = shuffleCijs(cijMeans,bestxMean,perturbType,thresh)
% [cijMeansShuff,cijStdsShuff] = shuffleCijs(cijMeans,bestxMean,perturbType,thresh)
% randomly perturb estimated cijMeans per state and calculate correct means
% perturbType='perState'
% perturbType='complete' % Might be very time consuming...
%
% Thomas Sauter, University of Luxembourg, 02/2011, thomas.sauter@uni.lu
% (c) 2014 University of Luxembourg Faculty of Science, Technology and Communication FSTC
% All rights reserved
% GPL version 3.0 to be found at: http://www.gnu.org/licenses/gpl.html

global estim
estim.display=0;
keep=estim.stopAtFirst;
estim.stopAtFirst=0;

orders=[];
diffs=[];
paramCounter=1;

switch perturbType
    case 'perState'
        NrRuns=sum(factorial(estim.paramNrList));
        runCounter=1;
        
        for counter=1:length(estim.optimizeStates)
            allPerms=perms((paramCounter+estim.paramNrList(counter)-1):-1:paramCounter)
            for counter2=1:size(allPerms,1)
                order=1:length(bestxMean);
                order(paramCounter:(paramCounter+estim.paramNrList(counter)-1))=allPerms(counter2,:)
                disp(['Cij perturbations: ' num2str(runCounter) ' / ' num2str(NrRuns)])
                diff=opt_fun(bestxMean(order));
                orders=[orders; order];
                diffs=[diffs; diff];
                runCounter=runCounter+1;
            end
            paramCounter=paramCounter+estim.paramNrList(counter);
        end
        
    case 'complete'
        NrRuns=factorial(sum(estim.paramNrList));
        runCounter=1;
        allPerms=perms(sum(estim.paramNrList):-1:1);
        for counter2=1:size(allPerms,1)
            order=allPerms(counter2,:)
            disp(['Cij perturbations: ' num2str(runCounter) ' / ' num2str(NrRuns)])
            diff=opt_fun(bestxMean(order));
            orders=[orders; order];
            diffs=[diffs; diff];
            runCounter=runCounter+1;
        end
        paramCounter=paramCounter+estim.paramNrList(counter);
end

disp('Done!')

estim.display=1;
estim.stopAtFirst=keep;

%
paramCounter=1;
cijsAllPert=[];
ok=find(diffs<thresh);
NrRunsIncluded=length(ok)
for counter=1:NrRunsIncluded %get all cijs of ok (below threshold)
    cijsPert=evalCijAndSimulate(estim,[],0,bestxMean(orders(ok(counter),:)));
    cijsAllPert(:,:,:,counter)=cijsPert; %not working for different sizes of cijs (per experiment)
    % solution: add additional rows with -1s
end
for counter=1:size(cijsAllPert,3) %display initial mean cij per experiment
    disp(['Starting guess EXPERIMENT' num2str(counter) ':'])
    disp('mean:')
    disp(mean(cijsAllPert(:,:,counter,:),4))
    disp('std:')
    disp(std(cijsAllPert(:,:,counter,:),0,4))
end
cijMeansPert=mean(cijsAllPert,4);
cijStdsPert=std(cijsAllPert,0,4);

cijMeansShuff=cijMeansPert;
cijStdsShuff=cijStdsPert;
paramCounter=1;
permsCounter=0;
for counter=1:length(estim.optimizeStates) %for all optimizeStates
    cijsAllPert=[];
    permsCounterOld=permsCounter;
    allPerms=perms((paramCounter+estim.paramNrList(counter)-1):-1:paramCounter);
    permsCounter=permsCounter+size(allPerms,1);
    ok_spec=ok(ok>permsCounterOld&ok<=permsCounter);
    disp(['Shuffling ' cell2mat(estim.optimizeStates(counter)) ' : '])
    NrRunsIncluded=length(ok_spec)
    for counter2=1:NrRunsIncluded %get state specific cijs (below threshold)
        cijsPert=evalCijAndSimulate(estim,[],0,bestxMean(orders(ok_spec(counter2),:)));
        cijsAllPert(:,:,:,counter2)=cijsPert;
    end
    
    cijMeansPert_spec=mean(cijsAllPert,4);
    cijStdsPert_spec=std(cijsAllPert,0,4);
    for counter2=1:estim.NrExps %check if optimizeState is optimized in specific experiment
        eval(['optimizeStatesExp=estim.exp' num2str(counter2) '.model.state_names(estim.exp' num2str(counter2) '.NrsOptimizeStates);']);
        eval(['state_namesExp=estim.exp' num2str(counter2) '.model.state_names;']);
        if ~isempty(intersect(estim.optimizeStates(counter),optimizeStatesExp))
            [C,IA,stateNr] = intersect(estim.optimizeStates(counter),state_namesExp);
            cijMeansShuff(:,stateNr,counter2)=cijMeansPert_spec(:,stateNr,counter2);
            cijStdsShuff(:,stateNr,counter2)=cijStdsPert_spec(:,stateNr,counter2);
        end
    end
    paramCounter=paramCounter+estim.paramNrList(counter);
end

disp('SHUFFLING RESULTS:')
for counter2=1:size(cijMeansShuff,3) %display
    disp(['EXPERIMENT' num2str(counter2) ':'])
    disp('mean:')
    disp(cijMeansShuff(:,:,counter2))
    %     disp('std:')
    %     disp(cijStdsShuff(:,:,counter2))
end

end

% --- End of script --- %