function [pdfs_ic,pdf_ic_mean,std_pdf_ic_mean,nsteps]=pbnStationary_TS(ics,n,F,varF,nf,nv,cij,state_names,eval_ss_params)
% [pdfs_ic,pdf_ic_mean,std_pdf_ic_mean,nsteps]=pbnStationary_TS(ics,n,F,varF,nf,nv,cij,state_names,eval_ss_params)
% Perform empirical simulation based calculation of the stationary distribution for a PBN
% Support 2 approaches of analysis
% 1) Two-state Markov chain approach
% - Check for convergence indicated by burn-in time (m0) and required time
% steps to reach steady-state distribution with defined confidence (N)
% 2) Random initial condition (multiple selections & fixed time steps)
% - randomly select initial conditions and calculate pdf (probability
% density function) based on the second half of Markov Chain
%
% Note: 3 functions from the BN/PBN toolbox are included in this function to increase efficiency (speed)
% "pbnRun.m", "pbnNextState.m", "pmfrnd.m", see Line 96 to 237
% Please refer to the scripts in the BN/PBN toolbox for additional information
%
% Thomas Sauter, University of Luxembourg, 09/2010, thomas.sauter@uni.lu
% Panuwat Trairatphisan, University of Luxembourg, 11/2014, panuwat.trairatphisan@uni.lu
% (c) 2014 University of Luxembourg Faculty of Science, Technology and Communication FSTC
% All rights reserved
% GPL version 3.0 to be found at: http://www.gnu.org/licenses/gpl.html

% Initialize parameters to store results
icsResult=[];
pdfs_ic=[];

% Retrieve the selected approach to analyze steady-state distribution
Approach=eval_ss_params(1);

if Approach==1 % Two-state Markov chain approach
    
    %  Retrive information on Two-state Markov Chain and perturbation parameters
    epsilon=eval_ss_params(2); % range of transition probability [Default=0.001]
    r=eval_ss_params(3); % range of accuracy (most sensitive) [Default=0.01 (IS) or =0.005 (RA & LS)]
    s=eval_ss_params(4); % probability of accuracy [Default=0.95]
    p=eval_ss_params(5); % perturbation in Shmulevich's scheme [Default=0, 0%]
    pMir=eval_ss_params(6); % perturbation in Miranda's & Parga's scheme [Deafault=0.001, 0.1%]
    initial_nsteps=eval_ss_params(7); % initial nsteps
    
    % Retrive information on output parameters
    optStates=eval_ss_params(8:end);
    
    % Initialize parameters
    m0=0; % initial burn-in steps
    nsteps=initial_nsteps; % initial MC simulation steps
    N_memory=0;
    m0_memory=0;
    N_collect=zeros(1,length(optStates));
    m0_collect=zeros(1,length(optStates));
    
    counting=1;
    y_kept=[];
    Collection=[];    
    pdf_ic=[];
    y=[];
    
    % Check for all inputs --> NO PERTURBATION!
    % Algorithm: if sum(nf*nv) >= 1 --> Indicator for input(s)
    % sum(nf*nv) = 0 --> input(s) with numerical values
    % sum(nf*nv) = 1 & varF = 'self-index' --> input which depends only on itself
    % sum(nf*nv) > 1 --> state with multiple functions --> NOT AN INPUT
    
    inputs_index=[];
    running_index=0;
    
    for counter=1:length(nf)
        nfnv=[];
        for counter2=1:nf(counter)
            nfnv=[nfnv; nf(counter)*nv(counter2+running_index)];
        end
        if sum(nfnv)==0 || (sum(nfnv)==1 && varF(1,counter)==counter)  
            inputs_index=[inputs_index counter];
        end
        running_index=running_index+nf(counter);
    end
    
    % Analyze for suitable m0 and N
            
    while size(y_kept,1)<nsteps
        
        % Simulation
        
        if counting==1 % First round of simulation
            
            Collection=[Collection; N_collect,m0_collect,N_memory,m0,nsteps];
            
            y=rand(1,n)>0.5; % random initial condition (Markov chain is already ergodic)

            if exist('ics')
                if ~isempty(ics) %fix specific initial conditions
                    for counter=1:size(ics,1)
                        y(ics(counter,1))=ics(counter,2);
                    end
                end
            end
            
            ynextrun=y; % Initialize next time step           
            
            y_keeping=zeros(nsteps,n); % pre-allocate matrix
            
            for counter=1:nsteps

                % Note: The functions called in the "quoted" command below are originally taken from the BN/PBN toolbox
                % All scripts are combined here to improve efficiency of the overall optPBN pipeline

                % "ynextrun = pbnRun(ynextrun(end,:),F,varF,nf,nv,cij,p,1);"
                
                % === pbnRun (from BN/PBN toolbox) === %
                % function y = pbnRun(x,F,varF,nf,nv,cij,p,nsteps)                
                
                x_pbnRun=ynextrun(end,:);
                nsteps_pbnRun=1;
                                
                n_pbnRun = length(nf); % number of genes
                
                if isstr(x_pbnRun) % if we want a random start
                    x_pbnRun = rand(1,n_pbnRun)>0.5;
                end                                
                
                y_pbnRun = zeros(nsteps_pbnRun+1,n_pbnRun); % initialize y
                y_pbnRun(1,:) = x_pbnRun; % first one is initial state
                
                for step=2:nsteps_pbnRun+1,

                    % "x = pbnNextState(x,F,varF,nf,nv,cij,p); %otherwise use logical functions"
                                        
                    % === pbnNextState (from BN/PBN toolbox) === %
                    % function y = pbnNextState(x,F,varF,nf,nv,cij,p)

                    x_pbnNextState=x_pbnRun;                    
                    
                    n_pbnNextState = length(x_pbnNextState);          % number of genes
                    y_pbnNextState = zeros(size(x_pbnNextState));     % initialize next state
                    gam = rand(1,n_pbnNextState) < p;    % generate a random perturbation vector
                    
                    if any(gam) % if gam is not all zeros
                        y_pbnNextState = bitxor(x_pbnNextState,gam); % perturb some genes
                    else
                        cnf = [0,cumsum(nf)];
                        b = 2.^[size(varF,1)-1:-1:0]';
                        
                        % Update each node separately.
                        for i=1:n_pbnNextState
                            pmf = cij(1:nf(i),i); % extract its cij probabilities
                                                        
                            % "j = pmfrnd(pmf); % pick a predictor at random (note: 1 <= j <= nf(i))"
 
                            % === pmfrnd (from BN/PBN toolbox) === % 
                            % function y = pmfrnd(pmf)
                            
                            cdf = cumsum(pmf);			% make cumulative distr. function
                            u = rand;						% pick a random number in [0,1]
                            j = sum(u > cdf) + 1;		% find where it falls
                                                        
                            k = cnf(i)+j; % Index of the random selected predictor.
                                                        
                            y_pbnNextState(i) = F(x_pbnNextState(varF(1:nv(k),k))*b(end-nv(k)+1:end)+1,k); % Output value for the selected function.
                        end % for i=1:n
                    end % if any(gam)
                    
                    ynextrun(step,:) = y_pbnNextState; % update history
                end
                                
                % === End of quoted scripts === %                
                
                pert=rand(1,n)<=pMir;
                pert(inputs_index)=0; % DO NOT PERTURB INPUTS
                ynext=mod(ynextrun(end,:)+pert,2);
                y_keeping(counter,:)=ynext;                            
            end
            
            y=y_keeping; % return trajectory
            
        elseif counting>1 % Next round of simulation
            
            y=y(end,:); % pick the last state from last simulation
            ynextrun=y; % initialize next state

            y_keeping=zeros(nsteps-size(y_kept,1)+1,n); % pre-allocate matrix            
            y_keeping(1,:)=y; % assign last state in the pre-allocated matrix
            
            for counter=size(y_kept,1):nsteps
                
                % Note: The functions called in the "quoted" command below are originally taken from the BN/PBN toolbox
                % All scripts are combined here to improve efficiency of the overall optPBN pipeline

                % "ynextrun = pbnRun(ynextrun(end,:),F,varF,nf,nv,cij,p,1);"
                
                % === pbnRun (from BN/PBN toolbox) === %
                % function y = pbnRun(x,F,varF,nf,nv,cij,p,nsteps)                
                
                x_pbnRun=ynextrun(end,:);
                nsteps_pbnRun=1;
                                
                n_pbnRun = length(nf); % number of genes
                
                if isstr(x_pbnRun) % if we want a random start
                    x_pbnRun = rand(1,n_pbnRun)>0.5;
                end                                
                
                y_pbnRun = zeros(nsteps_pbnRun+1,n_pbnRun); % initialize y
                y_pbnRun(1,:) = x_pbnRun; % first one is initial state
                
                for step=2:nsteps_pbnRun+1,

                    % "x = pbnNextState(x,F,varF,nf,nv,cij,p); %otherwise use logical functions"
                                        
                    % === pbnNextState (from BN/PBN toolbox) === %
                    % function y = pbnNextState(x,F,varF,nf,nv,cij,p)

                    x_pbnNextState=x_pbnRun;                    
                    
                    n_pbnNextState = length(x_pbnNextState);          % number of genes
                    y_pbnNextState = zeros(size(x_pbnNextState));     % initialize next state
                    gam = rand(1,n_pbnNextState) < p;    % generate a random perturbation vector
                    
                    if any(gam) % if gam is not all zeros
                        y_pbnNextState = bitxor(x_pbnNextState,gam); % perturb some genes
                    else
                        cnf = [0,cumsum(nf)];
                        b = 2.^[size(varF,1)-1:-1:0]';
                        
                        % Update each node separately.
                        for i=1:n_pbnNextState
                            pmf = cij(1:nf(i),i); % extract its cij probabilities
                                                        
                            % "j = pmfrnd(pmf); % pick a predictor at random (note: 1 <= j <= nf(i))"
 
                            % === pmfrnd (from BN/PBN toolbox) === % 
                            % function y = pmfrnd(pmf)
                            
                            cdf = cumsum(pmf);			% make cumulative distr. function
                            u = rand;						% pick a random number in [0,1]
                            j = sum(u > cdf) + 1;		% find where it falls
                                                        
                            k = cnf(i)+j; % Index of the random selected predictor.
                                                        
                            y_pbnNextState(i) = F(x_pbnNextState(varF(1:nv(k),k))*b(end-nv(k)+1:end)+1,k); % Output value for the selected function.
                        end % for i=1:n
                    end % if any(gam)
                    
                    ynextrun(step,:) = y_pbnNextState; % update history
                end
                                
                % === End of quoted scripts === %                

                pert=rand(1,n)<=pMir;
                pert(inputs_index)=0; % DO NOT PERTURB INPUTS
                ynext=mod(ynextrun(end,:)+pert,2);
                y_keeping(counter-size(y_kept,1)+2,:)=ynext;
            end

            y=y_keeping; % return trajectory
                        
        end
        
        % Append and keep the simulated trajectory 
        y_kept=[y_kept; y(2:end,:)]; % keep trajectory and discard IC/last state
        y=y_kept; % return trajectory for further calculation
        
        % Calculate nsteps (N) and m0 based on marginal distribution of each state
        [nsteps,m0,Collection] = TwoStateMC_MargDist(y,nsteps,m0,epsilon,r,s,optStates,Collection,N_memory,m0_memory,N_collect,m0_collect);
        
        % If the trajectory is not converged in 100 runs, something might go wrong
        counting=counting+1;
        
        if counting>100
            break
        end
        
    end
    
    % Reporting identified nsteps in the current iteration
    disp(['Number of nsteps = ' num2str(size(y,1)) ]);
    % Summary_Report=[Heading;num2cell(real(Collection))]

    % Calculate pdf for each state
    pdf_ic=mean(y(1:end,:),1);
    pdfs_ic=[pdfs_ic; pdf_ic];
        
elseif Approach==2
    
    p=0; % no perturbation
    
    %  Retrive information on Monte-Carlo simulation
    nsteps=eval_ss_params(2); % range of transition probability [Default=0.001]
    iterations=eval_ss_params(3); % range of accuracy (most sensitive) [Default=0.01 (IS) or =0.005 (RA & LS)]
    
    for counter=1:iterations
        x=rand(1,n)>0.5;
        if exist('ics')
            if ~isempty(ics) %fix specific initial conditions
                for counter=1:size(ics,1)
                    x(ics(counter,1))=ics(counter,2);
                end
            end
        end
        icsResult=[icsResult; x];
        pdf_ic=[];
        y = pbnRun(x,F,varF,nf,nv,cij,p,nsteps);
        pdf_ic=mean(y(floor(nsteps/2):end,:),1); % take 2nd half of simulation
        pdfs_ic=[pdfs_ic; pdf_ic];
    end
    
end

% Calculate mean and S.D. of pdf for each state
pdf_ic_mean=mean(pdfs_ic,1);
std_pdf_ic_mean=std(pdfs_ic,1);

end

% --- End of script --- %