function [X,FVAL] = genericEvA2SB(FUN,X,OPTIONS,method)
% genericEvA2SB: Generic interface to the EvA2 Interface, allowing for the
% following optimization algorithms (depending on the "method" setting (1-10):
%
% 1: Standard ES 
% 2: CMA-ES 
% 3: GA 
% 4: PSO 
% 5: DE 
% 6: Tribes 
% 7: Random (Monte Carlo) 
% 8: Hill-Climbing 
% 9: Cluster-based niching ES 
% 10: Clustering Hill-Climbing
%
% This function is not supposed to be called by itself, but rather from the
% optimization functions in the parent folder.
%
% Algorithms part of the EvA2 optimization package. Reference:
% Streichert F, Ulmer H. 2005. JavaEvA: A Java based framework for
% Evolutionary Algorithms, Tech. Report WSI-2005-06. Wilhelm-Schickard-
% Institut für Informatik (WSI), Centre for Bioinformatics Tübingen
% (ZBIT), Eberhard-Karls-University Tübingen, Germany.
%
% The interface from EvA2 to MATLAB was written by:
% Marcel Kronfeld
% Dept. of Computer Architecture
% University of Tübingen
% Germany 
%
% In case of a publication that was prepared using stdesJESB, users
% are asked to cite the above publication.
%
% NOTE: DO NOT BREAK THE OPTIMIZATION USING CTRL-C! If you do that the Java
% optimization thread will keep running and taking a lot of CPU performance
% away. Ultimately you will have to restart MATLAB. INSTEAD, use the
% button, which is provided by this function, to stop the optimization.
% In the case this function is used together with the parameter esimation
% functionality from the SBPD package the default button will be replaced
% by the SBPD parameter estimations own button.
%
% USAGE:
% ======
% [X,FVAL] = genericJESB(FUN,X,OPTIONS)
%
% FUN:      Function to optimize 
% X:        Starting Guess
% OPTIONS:  Structure containing options for the algorithm. 
%        OPTIONS.maxfunevals: Maximum number of function evaluations
%        OPTIONS.tolfun: Termination tolerance on the function value
%        OPTIONS.tolx: Termination tolerance on X
%        OPTIONS.highbounds: vector containing upper bounds for parameters
%           Instead of a vector highbounds can also be a scalar > 1. In the
%           latter case the highbounds are determined by multiplying this
%           scalar to the initial parameter guess X. This latter case works
%           only for positive X(i) ... if negative present then the user
%           needs to specify a whole highbounds vector. 
%        OPTIONS.lowbounds: vector containing lower bounds for parameters.
%           Instead of a vector lowbounds can also be a scalar < 1. In the
%           latter case the lowbounds are determined by multiplying this
%           scalar to the initial parameter guess X. This latter case works
%           only for positive X(i) ... if negative present then the user
%           needs to specify a whole lowbounds vector. 
%        OPTIONS.silent: =0: output of info, =1: no output
% method: integer defining the optimization method to be used (see above) 
%
% DEFAULT VALUES:
% ===============
% OPTIONS.maxfunevals:    2000*numberVariables
% OPTIONS.tolx:           0 (best to keep on 0!)
% OPTIONS.tolfun:         0 (best to keep on 0!)
% OPTIONS.lowbounds:      0.1  => lowbounds = 0.1*X 
% OPTIONS.highbounds:     10  => highbounds = 10*X 
% OPTIONS.silent:         0 (no output of info)
%
% Output Arguments:
% =================
% X:        Found solution (The best particle obtained for the population (Leader))
% FVAL:     Value of the function FUN at X 

% Information:
% ============
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  
% USA.

if method < 1 || method > 10,
    error('genericJESB: wrong argument value for ''method''.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONVERT FUNCTION NAME TO HANDLE (if necessary)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ischar(FUN),
    FUN = str2func(FUN);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ndim = length(X);
maxfunevals = 2000*ndim;
tolx = 0;
tolfun = 0;
lowbounds = 0.1*X;
highbounds = 10*X;
silent = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handle OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% silent
if isfield(OPTIONS,'silent'),
    if ~isempty(OPTIONS.silent),
        silent = OPTIONS.silent;
    end
end
% tolfun
if isfield(OPTIONS,'tolfun'),
    if ~isempty(OPTIONS.tolfun),
        tolfun = OPTIONS.tolfun;
    end
end
% tolx
if isfield(OPTIONS,'tolx'),
    if ~isempty(OPTIONS.tolx),
        tolx = OPTIONS.tolx;
    end
end
% maxfunevals
if isfield(OPTIONS,'maxfunevals'),
    if ~isempty(OPTIONS.maxfunevals),
        maxfunevals = OPTIONS.maxfunevals;
    end
end
% low and highbounds:
[lowbounds, highbounds] = handleLowHighBoundsSB(OPTIONS,X,lowbounds,highbounds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare and Call EvA2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('In case of a publication that was prepared using one of');
disp('the EvA2 optimizers, users are asked to cite the following publication:');
disp('Streichert F, Ulmer H. 2005. JavaEvA: A Java based framework for');
disp('Evolutionary Algorithms, Tech. Report WSI-2005-06. Wilhelm-Schickard-');
disp('Institut für Informatik (WSI), Centre for Bioinformatics Tübingen');
disp('(ZBIT), Eberhard-Karls-University Tübingen, Germany.');

% For some strange reason does the call to javaaddpath clear all global
% variables. Therefore we will store all globals locally first, add the
% classpath and finally write the globals back ... stupid but true!
allglobals = who('global');
for k=1:length(allglobals),
    eval(sprintf('global %s',allglobals{k}));
    eval(sprintf('global_save_%s = %s;',allglobals{k},allglobals{k}));
end

% ADD JE2Base.jar to the classpath
warning off;
javaaddpath(which('JE2Base.jar'));
warning on;

% Write back the global variables ... stupid stupid stupid!
for k=1:length(allglobals),
    eval(sprintf('global %s',allglobals{k}));
    eval(sprintf('%s = global_save_%s;',allglobals{k},allglobals{k}));
end

% Create optimset structure 
optionsJI = optimset('TolFun', tolfun);
optionsJI = optimset(optionsJI,'TolX', tolx);
optionsJI = optimset(optionsJI,'MaxFunEvals', maxfunevals);
if silent == 0,
    optionsJI = optimset(optionsJI,'Display', 'iter');
else
    optionsJI = optimset(optionsJI,'Display', 'off');
end

% Create JEInterface
R = [lowbounds; highbounds];
JI = JEInterface('JI',FUN,R,optionsJI); 

% Do the optimization
warning off;
JI = optimize(JI,method);
warning on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finish up and return the result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the result
[X,FVAL] = getResult(JI);    
X = X(:)'; % should be a row vector
return


