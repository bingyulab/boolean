function estim = add2estim(estim,n,nf,nv,F,varF,cij,state_names,optimizeStates,ics,measStates,meas,rules)
% estim = add2estim(estim,n,nf,nv,F,varF,cij,state_names,optimizeStates,ics,measStates,meas,rules)
% store PBN models with corresponding experimental data to estim global parameter
% generate the common list of optimize states and qualitative boundaries (if assigned)
%
% Thomas Sauter, University of Luxembourg, 02/2011, thomas.sauter@uni.lu
% Panuwat Trairatphisan, University of Luxembourg, 06/2014, panuwat.trairatphisan@uni.lu
% (c) 2014 University of Luxembourg Faculty of Science, Technology and Communication FSTC
% All rights reserved
% GPL version 3.0 to be found at: http://www.gnu.org/licenses/gpl.html

% Find the number of optimize states
[temp,NrsOptimizeStates,temp2]=intersect(state_names,optimizeStates);

for counter=length(optimizeStates):-1:1
    if nf(NrsOptimizeStates(counter))==1 %only 1 interaction function: not optimized
        NrsOptimizeStates(counter)=[];
    end
end

NrsOptimizeStates=sort(NrsOptimizeStates)';
optimizeStates=state_names(NrsOptimizeStates);

% Assign parameters to a single global parameter 'estim'
estim.NrExps=estim.NrExps+1;
temp=sum(sum(cij(:,NrsOptimizeStates)~=-1));
eval(['estim.exp' num2str(estim.NrExps) '.paramNr=temp;'])
eval(['estim.exp' num2str(estim.NrExps) '.NrsOptimizeStates=NrsOptimizeStates;'])
eval(['estim.exp' num2str(estim.NrExps) '.model.n=n;'])
eval(['estim.exp' num2str(estim.NrExps) '.model.nf=nf;'])
eval(['estim.exp' num2str(estim.NrExps) '.model.nv=nv;'])
eval(['estim.exp' num2str(estim.NrExps) '.model.F=F;'])
eval(['estim.exp' num2str(estim.NrExps) '.model.varF=varF;'])
eval(['estim.exp' num2str(estim.NrExps) '.model.cij=cij;'])
eval(['estim.exp' num2str(estim.NrExps) '.model.state_names=state_names;'])
eval(['estim.exp' num2str(estim.NrExps) '.meas.States=measStates;'])
eval(['estim.exp' num2str(estim.NrExps) '.meas.value=meas(1,:);'])
if size(meas,1)>=2
    eval(['estim.exp' num2str(estim.NrExps) '.SD.States=measStates;'])
    eval(['estim.exp' num2str(estim.NrExps) '.SD.value=meas(2,:);'])
end

% Collect information on qualitative bounds (if existed)
if size(rules,2)>2 % Once flags are assigned
    UB=[];LB=[];
    for counter=1:size(rules,1)
        
        if strcmp(cell2mat(rules(counter,3)),'C') % Constant/Fixed
            UB=UB; LB=LB;            
        
        elseif strcmp(cell2mat(rules(counter,3)),'D') % Default/Dynamic/Free
            UB=[UB 1]; LB=[LB 0];            
        
        elseif strcmp(cell2mat(rules(counter,3)),'H') % Higher weight
            UB=[UB 1]; LB=[LB estim.BoundCutoff];            

        elseif strcmp(cell2mat(rules(counter,3)),'L') % Lower weight
            UB=[UB estim.BoundCutoff]; LB=[LB 0];
            
        end
    end
    eval(['estim.exp' num2str(estim.NrExps) '.UB=UB;'])
    eval(['estim.exp' num2str(estim.NrExps) '.LB=LB;'])
end
eval(['estim.exp' num2str(estim.NrExps) '.ics=ics;'])
eval(['estim.exp' num2str(estim.NrExps) '.rules=rules;'])

end

% --- End of script --- %