function res = binNum2dec(y,n)
% res = binNum2dec(y,n)
% convert binary vector y to decimal + 1!!! [for BNPBN toolbox]
% n: number of bits
%
% Thomas Sauter, University of Luxembourg, 09/2010, thomas.sauter@uni.lu
% (c) 2014 University of Luxembourg Faculty of Science, Technology and Communication FSTC 
% All rights reserved
% GPL version 3.0 to be found at: http://www.gnu.org/licenses/gpl.html

b = 2.^[n-1:-1:0]'; % Used for binary / decimal conversion

res = y*b + 1; %convert vector y to bin

end

% --- End of script --- %