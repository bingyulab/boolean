function diff = opt_fun(x)
% diff = opt_fun(x)
% evaluate optimal cost/objective function from the set of sampling parameters via simulation
% optimal cost = sum of square error (SSE) or chi-square between predicted state values from simulation and measurement data
% results (optimal cost and respective parameters) are saved in separate file (tempFIT)
%
% Thomas Sauter, University of Luxembourg, 02/2011, thomas.sauter@uni.lu
% Panuwat Trairatphisan, University of Luxembourg, 06/2014, panuwat.trairatphisan@uni.lu
% (c) 2014 University of Luxembourg Faculty of Science, Technology and Communication FSTC
% All rights reserved
% GPL version 3.0 to be found at: http://www.gnu.org/licenses/gpl.html

global estim

% Retrieve parameters from estim
paramNr=estim.paramNr;
paramNrList=estim.paramNrList;
allParamNrList=estim.paramNrList;
optimizeStates=estim.optimizeStates;
fileName=estim.fileName;
level=estim.optLevels;
Approach=estim.Approach;
if isfield(estim,'objFun') == 1
    objFun=estim.objFun;
end

% Initialize parameters to calculate SSE and to keep results
diff=0;
history_new=[];

% Extract information from estim and evaluate functions
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
    
    % Generate common list of optimize states
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
        
    % Calculate difference between model simulation and measurement data
    if isfield(estim,'objFun') ~= 1
        diff=diff+sum((pdf_ic_mean(NrsMeasStates)-meas(orderMeas)).^2);
    else    
        if objFun==1
            diff=diff+sum((pdf_ic_mean(NrsMeasStates)-meas(orderMeas)).^2);
        elseif objFun==2 && ~isempty(SD)
            diff=diff+sum(((pdf_ic_mean(NrsMeasStates)-meas(orderMeas)).^2)./SD(orderMeas).^2);
        elseif objFun==2 && isempty(SD)
            warning('Warning: Please provide SD value to calculate chi-square value')
            pause(1)
        else
            error('Error: Please review the choice of your objective function (SSE or chi-square)')
        end
    end
        
    % Display intermediate results during the optimization
    if estim.display
        cij
        [state_names(NrsMeasStates) num2cell(pdf_ic_mean(NrsMeasStates)')]
        diff
    end
    history_new=[history_new, pdf_ic_mean(NrsMeasStates)];
    
end

% Store results into a file
history_new=[diff, history_new, x];

if exist(fileName)
    load(fileName,'history','history_old')
else
    history=[]; history_old=[];
end

history=[history; history_new];
if size(history,1)>9999
    history_old=history;
    history=[];
end
save(fileName)

% If good quality fit is reached, stop optimization
if (diff<=estim.thresh) && estim.stopAtFirst
    error('First good fit reached')
end
% toc
end

% --- End of script --- %