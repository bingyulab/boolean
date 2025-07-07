function [n,nf,nv,F,varF,cij,state_names,UB,LB] = rule2PBN(rules)
% [n,nf,nv,F,varF,cij,state_names] = rule2PBN(rules)
% construct PBN model from logical rule list
%
% Examples:
% rules={};
% rules=[rules; {'TNFa = TNFa','1'}];
% rules=[rules; {'VD3 = VD3 | ~TNFa','1'}];
% rules=[rules; {'VDR = VD3 & TNFa','0.4'}];
% rules=[rules; {'VDR = VD3','0.6'}];
%
% numerical rules:
% rules=[rules; {'N1 = 1','1'}];
% rules=[rules; {'N1 = 1','0.8'}]; %no need to specify 'N1 = 0','0.2' in addition!!!
% rules=[rules; {'N1 = 0','1'}]; %no concentration inputs for '=0' !!!
%
% combined variables and numercal rules (ordering: concentration inputs last!):
% rules=[rules; {'N1 = N1','0.5'}];
% rules=[rules; {'N1 = 1','0.5'}];
% or:
% rules=[rules; {'N1 = N1','0.5'}];
% rules=[rules; {'N1 = 0','0.5'}];
%
% Thomas Sauter, University of Luxembourg, 09/2010, thomas.sauter@uni.lu
% Panuwat Trairatphisan, University of Luxemoburg, 06/2014, panuwat.trairatphisan@uni.lu
% (c) 2014 University of Luxembourg Faculty of Science, Technology and Communication FSTC
% All rights reserved
% GPL version 3.0 to be found at: http://www.gnu.org/licenses/gpl.html

%identify all variables and inputs / outputs and number of variables in each Boolean function (nv)
variables=[];
inputs=[];
outputs=[];
nv_inputorder=[];
for counter=1:size(rules,1)
    rule=cell2mat(rules(counter));
    newVariables=regexp(rule,'\w*','match');
    newVariables=newVariables(find(~strcmp(newVariables,'1'))); %remove '1' from variable list
    newVariables=newVariables(find(~strcmp(newVariables,'0'))); %remove '0' from variable list
    variables=[variables, newVariables];
    outputs=[outputs, newVariables(1)];
    inputs_all=newVariables(2:end);
    [temp,I,temp2] =unique(inputs_all,'first');
    [temp,II]=sort(I);
    inputs_unique=inputs_all(I(II));
    inputs=[inputs, inputs_unique];
    
    %   works with Matlab7.7 (not 7.8)
    if ~isempty(newVariables(2:end))
        nv_inputorder=[nv_inputorder, size(unique(newVariables(2:end)),2)];
    else
        nv_inputorder=[nv_inputorder, 0];
    end
    if length(newVariables)==1 %numerical input
        conc=str2num(cell2mat(rules(counter,2)));
        if ~(conc==round(conc)) %non-integer numerical input ("concentration")
            nv_inputorder=[nv_inputorder, 0]; %add additional function
            outputs=[outputs, newVariables(1)];
        end
    end
end

[temp,I,temp2] =unique(variables,'first');
[temp,II]=sort(I);
variables=variables(I(II)); %!!!!!!!!!!!!!!!! not used: replaced by output_unique
% might cause problems if variables ~= output_unique (e.g. additional variables)

nv_inputorder;
[temp,I,temp2] =unique(outputs,'first');
[temp,II]=sort(I);
outputs_unique=outputs(I(II));
state_names=outputs_unique';
n=length(outputs_unique);
nv=[];
outputs_ordered=[];
outputs_ordering=[];
for counter=1:length(outputs_unique)
    for counter2=1:length(outputs)
        if (isequal(outputs(counter2),outputs_unique(counter)))
            nv=[nv, nv_inputorder(counter2)];
            outputs_ordered=[outputs_ordered, outputs(counter2)];
            outputs_ordering=[outputs_ordering, counter2];
        end
    end
end
nv;
outputs_ordered;

%nf =  number of functions per node
nf=zeros(length(outputs_unique),1)';
for counter=1:length(outputs_unique)
    for counter2=1:length(outputs)
        nf(counter)=nf(counter)+(isequal(outputs(counter2),outputs_unique(counter)));
    end
end
nf;

%cij = probabilities (rows) per output node (columns)
cij=-ones(max(nf),length(outputs_unique));
concInputCounter=0;
for counter=1:length(outputs_unique)
    counter3=0;
    counter2=1;
    while counter2<=length(outputs)
        if (isequal(outputs(counter2),outputs_unique(counter)))
            counter3=counter3+1;
            conc=str2num(cell2mat(rules(counter2-concInputCounter,2)));
            if conc==0
                cij(counter3,counter)=1;
            else
                cij(counter3,counter)=conc;
            end
            if nv(counter2)==0 %numerical input
                if ~(conc==round(conc)) %non-integer numerical input ("concentration")
                    concInputCounter=concInputCounter+1;
                    counter3=counter3+1;
                    cij(counter3,counter)=1-conc;
                    counter2=counter2+1;
                end
            end
        end
        counter2=counter2+1;
    end
end
cij;

%input variables
varF=-ones(max(nv),length(outputs));
varF(1,1:length(outputs))=ones(length(outputs),1); %input node1 for all numerical inputs (doesn't matter)
for counter=1:length(outputs)
    counter3=sum(nv_inputorder(1:(outputs_ordering(counter)-1)))+1;
    for counter2=1:nv(counter)
        [temp,IA,temp2] = intersect(outputs_unique,inputs(counter3));
        %         [temp,IA,temp2] = intersect(variables,inputs(counter3));
        varF(counter2,counter)=IA;
        counter3=counter3+1;
    end
end
if size(varF,1)==1 %!!!!!!!!! max nv=1 => pbnNextState does not work ("Panuwat, 7.2.11")
    varF(2,:)=-ones(1,sum(nf)); %add one more line
end
varF;
% variables
outputs_unique;

%interaction functions
F=-ones(2^max(nv),length(outputs));
concInputCounter=0;
counter=1;
while counter<=length(outputs) %for all outputs (incl +1/concentration)
    if nv(counter)~=0 %logical input
        for counter2=0:((2^nv(counter))-1) %length of interaction table
            BoolInput=dec2bin(counter2,nv(counter));
            for counter3=1:length(BoolInput) %walk through all possible input values
                eval([cell2mat(outputs_unique(varF(counter3,counter))) '=' num2str(BoolInput(counter3)) ';'])
            end
            rule=cell2mat(rules(outputs_ordering(counter)-concInputCounter)); %get active rule
            inputRuleStart=findstr(rule,'=')+1;
            inputRule=rule(inputRuleStart:end);
            %         eval(inputRule)
            F(counter2+1,counter)=eval(inputRule); %evaluate rule for specific input values
        end
    else %numerical input
        rule=cell2mat(rules(outputs_ordering(counter)-concInputCounter)); %get active rule
        inputRuleStart=findstr(rule,'=')+1;
        inputRule=rule(inputRuleStart:end);
        conc=str2num(cell2mat(rules(outputs_ordering(counter)-concInputCounter,2)));
        if (conc==0 | (eval(inputRule)==0))
            F(1:2,counter)=0; %always 0
        else
            F(1:2,counter)=1; %always 1
        end
        if ~(conc==round(conc)) %non-integer numerical input ("concentration")
            counter=counter+1;
            F(1:2,counter)=0;
            concInputCounter=concInputCounter+1;
        end
    end
    counter=counter+1;
end
F;

% allow for the combination of variable and concentration input
% => remove added 0 functions in this case
for counter=1:length(outputs_unique)
    if nf(counter)>2 %more than 2 functions per output
        specNrInputs=nv((sum(nf(1:counter-1))+1):sum(nf(1:counter)));
        if sum(specNrInputs==0)>1 %more than 1 numerical input
            %            "0"-functions added wrongly
            WhereZeros=find(specNrInputs==0);
            removeZero=sum(nf(1:counter-1))+WhereZeros(2); %identify index of zero which as to be removed
            
            nv(removeZero)=[]; %remove all unnecessary entries
            nf(counter)=nf(counter)-1;
            F(:,removeZero)=[];
            varF(:,removeZero)=[];
            cij(WhereZeros(2),counter)=-1;
            if max(cij(end,:))==-1
                cij(end,:)=[];
            end
        end
    end
end

% --- End of script --- %