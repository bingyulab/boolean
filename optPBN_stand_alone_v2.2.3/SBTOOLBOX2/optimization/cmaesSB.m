function [X,FVAL] = cmaesSB(varargin)
% cmaesSB: (15,50)-Evolution Strategy with Covariance Matrix Adaptation.
% (2nd of the EvA2 optimizers)
%
% Algorithm part of the EvA2 optimization package. Reference:
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
% In case of a publication that was prepared using cmaesSB, users
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
% [info] = cmaesSB()
% [X,FVAL] = cmaesSB(FUN,X)
% [X,FVAL] = cmaesSB(FUN,X,OPTIONS)
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
% info: calling the function w/o input argument returns information about
%       the options and a flag indicating if the algorithm can handle
%       constraints or not
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0,
    X = [];
    X.name = 'cmaesSB';
    X.constrained = 1;
    X.description = '(15,50)-Evolution Strategy with Covariance Matrix Adaptation (global)';          
    X.defaultOptions.names = {'maxfunevals', 'tolfun', 'tolx'};
    X.defaultOptions.values = {'50000','0','0'};
    X.defaultOptions.description = {'Maximum number of function evaluations', 'Termination tolerance on the function value', 'Termination tolerance on X'};
    FVAL = [];
    return
elseif nargin == 2,
    FUN = varargin{1};
    X = varargin{2};
    OPTIONS = [];
elseif nargin == 3,
    FUN = varargin{1};
    X = varargin{2};
    OPTIONS = varargin{3};
else
    error('Incorrect number of input arguments.');
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the generic EvA2 interface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
method = 2;
[X,FVAL] = genericEvA2SB(FUN,X,OPTIONS,method);    
X = X(:)'; % should be a row vector
return


