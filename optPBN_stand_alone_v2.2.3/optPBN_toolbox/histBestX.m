function histBestX(bestx,NrX)
% histBestX(bestx,NrX)
% generate histogram of NrX diff values(optimal cost/objective function) of bestx
%
% Thomas Sauter, University of Luxembourg, 02/2011, thomas.sauter@uni.lu
% (c) 2014 University of Luxembourg Faculty of Science, Technology and Communication FSTC
% All rights reserved
% GPL version 3.0 to be found at: http://www.gnu.org/licenses/gpl.html

diffs=[];
for counter=1:NrX
    diff=opt_fun(bestx)
    diffs=[diffs, diff];
end
figure
hist(diffs)

end

% --- End of script --- %