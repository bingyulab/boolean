function res = dec2binNum(y,n)
% res = dec2binNum(y,n)
% convert dec -1 !!! to binary vector [for BNPBN toolbox]
% input y: column vector of decimal numbers
% n: number of bits
%
% Thomas Sauter, University of Luxembourg, 09/2010, thomas.sauter@uni.lu
% (c) 2014 University of Luxembourg Faculty of Science, Technology and Communication FSTC
% All rights reserved
% GPL version 3.0 to be found at: http://www.gnu.org/licenses/gpl.html

bits = [n:-1:1];

res=[];
for counter=1:size(y,1)
    res = [res; bitget(y(counter)-1,bits)]; %convert number y to dec
end

end

% --- End of script --- %