function cijs = evalCijAndSimulate_G5K(estim,bestRun,simFlag,Memory,x)
% cijs = evalCijAndSimulate_G5K(estim,bestRun,simFlag,Memory,x)
% convert parameters to selection probability (cij) and simulate the model based on calculated Cij
% can use either bestRun number or opimized parameters x (bestRun=[])
% predicted state values from model simulation are compared to measurement data
% Note: The results from the simulations might not always be the same but are closely clustered
% Note: This script is customized for the post-processing of the result from grid optimization
% Grid5000/Cluster-based variance of opt_fun script
%
% Thomas Sauter, University of Luxembourg, 02/2011, thomas.sauter@uni.lu
% Panuwat Trairatphisan, University of Luxembourg, 06/2014, panuwat.trairatphisan.@uni.lu
% (c) 2014 University of Luxembourg Faculty of Science, Technology and Communication FSTC
% All rights reserved
% GPL version 3.0 to be found at: http://www.gnu.org/licenses/gpl.html

% Extract information from the global parameter estim
Approach=estim.Approach;
nsteps=estim.nsteps;
iterations=estim.iterations;
level=estim.optLevels;
paramNr=estim.paramNr;
allParamNrList=estim.paramNrList;
optimizeStates=estim.optimizeStates;
if isfield(estim,'objFun') == 1
    objFun=estim.objFun;
end

% Load results from file
if ~isempty(bestRun)
    x=Memory(bestRun, 3:paramNr+2);
end

acc_diff=0; %introduce accumulated differences (optimal cost/objective function)

cijs=[];
for counter3=1:estim.NrExps
    eval(['n=estim.exp' num2str(counter3) '.model.n;'])
    eval(['nf=estim.exp' num2str(counter3) '.model.nf;'])
    eval(['nv=estim.exp' num2str(counter3) '.model.nv;'])
    eval(['F=estim.exp' num2str(counter3) '.model.F;'])
    eval(['varF=estim.exp' num2str(counter3) '.model.varF;'])
    eval(['cij=estim.exp' num2str(counter3) '.model.cij;'])
    eval(['state_names=estim.exp' num2str(counter3) '.model.state_names;'])
    eval(['measStates=estim.exp' num2str(counter3) '.meas.States;'])
    eval(['meas=estim.exp' num2str(counter3) '.meas.value;'])
    if isfield(eval(['estim.exp' num2str(counter3)]),'SD') == 1
        eval(['SD=estim.exp' num2str(counter3) '.SD.value;'])
    else
        SD=[];    
    end
    if isfield(eval(['estim.exp' num2str(counter3)]),'UB') == 1 && isfield(eval(['estim.exp' num2str(counter3)]),'LB') == 1        
        eval(['UB=estim.exp' num2str(counter3) '.UB;'])
        eval(['LB=estim.exp' num2str(counter3) '.LB;'])
    end       
    eval(['paramNrExp=estim.exp' num2str(counter3) '.paramNr;'])
    eval(['NrsOptimizeStatesExp=estim.exp' num2str(counter3) '.NrsOptimizeStates;'])
    eval(['ics=estim.exp' num2str(counter3) '.ics;'])
    
    [temp,NrsMeasStates,orderMeas]=intersect(state_names,measStates);
    [StatesToOptimizeExp,NrsOptimizeStatesExp_orderedAll,temp2]=intersect(optimizeStates,state_names(NrsOptimizeStatesExp)');
    [temp,temp2,NrsOptimizeStatesExp]=intersect(StatesToOptimizeExp,state_names');
    
    paramNrCounter=0;
    paramNrCounterIncluded=0;
    
    % Assign sampling parameter (x) to the appropriate position of selection probability (cij) in the model
    for counter=1:length(optimizeStates)
        
        if max(counter==NrsOptimizeStatesExp_orderedAll) %state to be optimized in this particular experiment?
            paramNrCounterIncluded=paramNrCounterIncluded+1;
            if cij(end,NrsOptimizeStatesExp(paramNrCounterIncluded))==(-1)
                largestEntry=min(find(cij(:,NrsOptimizeStatesExp(paramNrCounterIncluded))==-1)-1);
            else
                largestEntry=size(cij,1);
            end
            for counter2=1:largestEntry %(largestEntry-1)
                paramNrCounter=paramNrCounter+1;

                if exist('UB') && exist('LB') % If bounds existed -> assign sampling parameter (x) based on bounds
                    pre_cij=x(paramNrCounter);
                    current_cij=(pre_cij*(UB(paramNrCounter)-LB(paramNrCounter))+LB(paramNrCounter));
                    cij(counter2,NrsOptimizeStatesExp(paramNrCounterIncluded))=current_cij;                    
                else % if bounds doesn't exist -> use default bounds ([0 1]) for all optimizing parameters                    
                    cij(counter2,NrsOptimizeStatesExp(paramNrCounterIncluded))=x(paramNrCounter);
                end
                
            end
        else
            paramNrCounter=paramNrCounter+allParamNrList(counter); %skip parameters per this state (not optimized here)
        end
    end
      
    % If discrete mode is chosen, discretize cij to discrete values
    if level>=1
        cij=max(-1,min(level-1,floor(cij*level))/(level-1)); % discritize to level
    end  
    
    % Normalize selection probability (cij) to the sum of all parameter values for each node -> final cij
    for counter=1:size(cij,2)
        if cij(end,counter)==(-1)
            largestEntry=min(find(cij(:,counter)==-1)-1);
        else
            largestEntry=size(cij,1);
        end
        if sum(cij(1:largestEntry,counter))==0 %only zero probabilities?
            cij(1:largestEntry,counter)=ones(largestEntry,1); %use all then
        end
        cij(1:largestEntry,counter)=cij(1:largestEntry,counter)./sum(cij(1:largestEntry,counter));
    end
    
    % Minor fix for inconsistent cijs
    if size(cij,1)~=max(size(cijs,1))
        ToAdd=ones((max(size(cijs,1))-size(cij,1)),size(cij,2)).*-1;
        cij=[cij; ToAdd];
    end
    
    cijs(:,:,counter3)=cij;

    if simFlag
        
        % Collect and prepare parameters for steady-state calculation
        optStates_index=(find(ismember(state_names,optimizeStates)))';
        if Approach==1
            epsilon=estim.epsilon;
            r=estim.r;
            s=estim.s;
            p=estim.p;
            pMir=estim.pMir;
            initial_nsteps=estim.nsteps;
            eval_ss_params=[Approach, epsilon, r, s, p, pMir, initial_nsteps, optStates_index];
        elseif Approach==2
            nsteps=estim.nsteps;
            iterations=estim.iterations;
            eval_ss_params=[Approach, nsteps, iterations];
        end
        
        % Calculate steady-state distribution
        [pdfs_ic,pdf_ic_mean,std_pdf_ic_mean,nsteps]=pbnStationary_TS(ics,n,F,varF,nf,nv,cij,state_names,eval_ss_params);
        measDisplay(1:length(state_names))=NaN; measDisplay(NrsMeasStates)=meas(orderMeas);
        disp(['EXPERIMENT ' num2str(counter3) ':']), disp(' ')
        disp([{'State','Model','Meas'}; state_names num2cell(pdf_ic_mean') num2cell(measDisplay')])
        
        % Generate the summation of diff value (aka optimal cost) from each experimental case
        if isfield(estim,'objFun') ~= 1
            acc_diff=acc_diff+sum((pdf_ic_mean(NrsMeasStates)-meas(orderMeas)).^2)
        else
            if objFun==1
                acc_diff=acc_diff+sum((pdf_ic_mean(NrsMeasStates)-meas(orderMeas)).^2)
            elseif objFun==2 && ~isempty(SD)
                acc_diff=acc_diff+sum(((pdf_ic_mean(NrsMeasStates)-meas(orderMeas)).^2)./SD(orderMeas).^2)
            elseif objFun==2 && isempty(SD)
                warning('Warning: Please provide SD value to calculate chi-square value')
                pause(1)
            else
                error('Error: Please review the choice of your objective function (SSE or chi-square)')
            end
        end

        
    end
end

end

% --- End of script --- %