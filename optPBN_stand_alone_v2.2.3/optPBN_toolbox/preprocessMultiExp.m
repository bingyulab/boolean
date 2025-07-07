function estim = preprocessMultiExp(estim)
% estim = preprocessMultiExp(estim)
% Preprocessing for multiple experiments
% combine parameters (and boundaries) of all experiments to a unique list
% & get Nr of necessary parameters for optimization process
%
% Thomas Sauter, University of Luxembourg, 02/2011, thomas.sauter@uni.lu
% Panuwat Trairatphisan, University of Luxembourg, 06/2014, panuwat.trairatphisan@uni.lu
% (c) 2014 University of Luxembourg Faculty of Science, Technology and Communication FSTC
% All rights reserved
% GPL version 3.0 to be found at: http://www.gnu.org/licenses/gpl.html

allParamsNUnique={};
allParamNrs=[];
for counter=1:estim.NrExps
    eval(['allParamsNUnique=[allParamsNUnique; estim.exp' num2str(counter) '.model.state_names(estim.exp' num2str(counter) '.NrsOptimizeStates)];']);
    eval(['allParamNrs=[allParamNrs, sum(estim.exp' num2str(counter) '.model.cij(:,estim.exp' num2str(counter) '.NrsOptimizeStates)~=-1)];']);
    sum(estim.exp2.model.cij(:,estim.exp1.NrsOptimizeStates)~=-1)
    allParamNrs(find(allParamNrs==0))=[];        
end

[allParams,I,J]=unique(allParamsNUnique);
estim.optimizeStates=allParams';

for counter=1:max(J) %length(allParams)
    if min(allParamNrs(J==counter))~=max(allParamNrs(J==counter))
        disp(['ERROR: Nr of interaction functions for node ' cell2mat(unique(allParamsNUnique(J==counter))) ' not matching for different experiments'])
        error('This is currently not supported by this Toolbox. Local variables necessary. [TS, 15.12.10]')
    end
end
estim.paramNr=sum(allParamNrs(I));
estim.paramNrList=allParamNrs(I); %parameters per optimizeStates

if isfield(estim.exp1,'UB') == 1 && isfield(estim.exp1,'LB') == 1 % if flags for qualitative boundaries are assigned
    
    % Ordering and sorting of bounds
    preOrder=sortrows([I estim.paramNrList']); % Order of optParam + Nr of each optParam (sampled parameter-based)
    
    for counter=1:size(preOrder,1)
        if counter==1
            counting=1;
            counting2=preOrder(1,2);
        else
            counting=counting+preOrder(counter-1,2); % StartPos: Adding up ParamNr of the previous row
            counting2=counting2+preOrder(counter,2); % EndPos
        end
        preOrder(counter,3)=counting;
        preOrder(counter,4)=counting2;
    end
    
    Order=zeros(size(preOrder,1),size(preOrder,2)); % Shuffle rows based on the order of optParams
    for counter=1:size(preOrder,1)
        Order(counter,:)=preOrder(I(counter),:);
    end
    
    % Assign correct boundaries to estim
    All_UB=[];All_LB=[]; % Initialize parameters to collect upper and lower qualitative bounds
    for counter=1:estim.NrExps
        % Extract upper and lower bounds from rules
        eval(['current_UB=[estim.exp' num2str(counter) '.UB];']);
        eval(['current_LB=[estim.exp' num2str(counter) '.LB];']);
        UB=zeros(1,length(current_UB));
        LB=zeros(1,length(current_LB));
        Bounds_counter=1;
        for counter2=1:size(Order,1) % for each optParam
            for counter3=Order(counter2,3):Order(counter2,4) % Search for StartPos:EndPos of sampled parameters
                UB(Bounds_counter)=current_UB(counter3); % Assign to UB
                LB(Bounds_counter)=current_LB(counter3); % Assign to LB
                Bounds_counter=Bounds_counter+1;
            end
        end
        eval(['[estim.exp' num2str(counter) '.UB]=UB;']);
        eval(['[estim.exp' num2str(counter) '.LB]=LB;']);
        % Check if the number of qualitative boundaries is consistent in all experiments
        try
            All_UB=[All_UB; UB];All_LB=[All_LB; LB]; % Appending bounds from each experiment
        catch
            disp('ERROR: The number of qualitative boundaries is not consistest across experiments')
            error(['Please verify the number of constant boundaries (C) in experiment (case) Nr: ' num2str(counter) ])
        end
    end

    % Check if the boundaries are consistent across different experiments
    % if size(unique(All_UB,'rows'),1) > 1 || size(unique(All_LB,'rows'),1) > 1  % Check all rows at once
    for counter=1:estim.NrExps-1 % Check current row with the next row 
        if  size(unique(All_UB(counter:counter+1,:),'rows'),1) > 1 || size(unique(All_LB(counter:counter+1,:),'rows'),1) > 1
            disp('ERROR: The assignment of qualitative boundaries are not consistest across experiments')
            error(['Please verify the choice of boundary setting in experiment (case) Nr: ' num2str(counter+1)])
        end
    end
end

disp(estim) % display final estim if there is no error

end

% --- End of script --- %