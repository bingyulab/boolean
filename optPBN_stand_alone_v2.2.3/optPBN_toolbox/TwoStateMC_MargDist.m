function [nsteps,m0,Collection] = TwoStateMC_MargDist(y,nsteps,m0,epsilon,r,s,optStates,Collection,N_memory,m0_memory,N_collect,m0_collect)
% [nsteps,m0,Collection] = TwoStateMC_MargDist(y,nsteps,m0,epsilon,r,s,optStates,Collection,N_memory,m0_memory,N_collect,m0_collect)
% Identify burn-in steps (m0) and required simulation steps with confience (N)
% Based on the evaluation of marginal distribution of each optimized state
%
% Thomas Sauter, University of Luxembourg, 01/2014, thomas.sauter@uni.lu
% Panuwat Trairatphisan, University of Luxembourg, 04/2014, panuwat.trairatphisan@uni.lu
% (c) 2014 University of Luxembourg Faculty of Science, Technology and Communication FSTC 
% All rights reserved
% GPL version 3.0 to be found at: http://www.gnu.org/licenses/gpl.html

% Calculate N and m0 from marginal distribution of each measured states

for counter=1:length(optStates)
    
    y_current=y(:,optStates(counter)); % select one state at a time
    
    % Can also use crosstab function to count states but it is slower
    % counts=crosstab(y_current(m0+1:end-1),y_current(m0+2:end));    

    y_crosstab=[y_current(m0+1:end-1),y_current(m0+2:end)];
    pre1=sum(y_crosstab(:,1));
    counts(2,2)=sum(y_crosstab(:,1).*y_crosstab(:,2));
    counts(2,1)=pre1-counts(2,2);
    counts(1,1)=sum((y_crosstab(:,1)-1).*(y_crosstab(:,2)-1));
    counts(1,2)=size(y_crosstab,1)-counts(1,1)-counts(2,1)-counts(2,2);
   
    pp=counts./repmat(sum(counts,2),1,2);
    % pp(1,2)/(pp(1,2)+pp(2,1)); % Quality control to check if the sum of each row is 1
       
    % Extract alpha (probability of leaving current state) and beta (probability of coming to current state)
    alpha=pp(1,2);
    beta=pp(2,1);
        
    % **Calculate m0 and N for each gene and store it
    m0_temp=log10(epsilon*(alpha+beta)/(max(alpha,beta)))/(log10(abs(1-alpha-beta)));
    N_temp=alpha*beta*(2-alpha-beta)/(alpha+beta)^3 * (r/norminv(0.5*(1+s)))^(-2);    
    % Note: The formulas to calculate m0 and N have been corrected.
    
    m0_collect(1,counter)=real(m0_temp);
    N_collect(1,counter)=N_temp;
    
end

% Check for the highest number of m0 and N and keep in memory
N_memory=ceil(max([max(N_collect),N_memory]));
m0_memory=ceil(max([max(m0_collect),m0_memory]));

% Only if more nsteps is required, nsteps and m0 will be updated
if nsteps<N_memory
    nsteps=N_memory;
    m0=m0_memory;
end

% Collect results for further reporting system
Collection=[Collection; N_collect,m0_collect,N_memory,m0,nsteps];

end

% --- End of script --- %