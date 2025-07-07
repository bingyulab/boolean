function [ab,d,attractor_means,attracProb_ic,pdf_ic,pmf_ic]=calcAttractors_ic(n,F,varF,nv,ics,state_names)
% function calcAttractors_ic(n,F,varF,nv,ic,state_names)
% calculate attractors and basins for specific initial conditions
%
% Thomas Sauter, University of Luxembourg, 09/2010, thomas.sauter@uni.lu

[A,Avec] = bnAsparse(F,varF,nv);
[ab,d] = bnAttractor(Avec);
s = bnStats(ab,d); %statistics of a Boolean Network
all_states_bin=dec2binNum([1:size(A,1)]',n);

ic_sel=ones(1,size(all_states_bin,1));
if exist('ics')
    if ~isempty(ics)
        for counter=1:size(all_states_bin,1)
            if ~(min(all_states_bin(counter,ics(:,1))==ics(:,2)'))
                ic_sel(counter)=0;
            end
        end
    end
end
disp([num2str(sum(ic_sel==1)) ' / ' num2str(length(ic_sel)) ' states selected as initial conditions'])

% display attractors etc.
summary=[];
attractor_means=[];
for counter=-1:-1:min(ab)
    disp(['attractor ' num2str(abs(counter)) ', basin size = ' num2str(sum(ab(find(ic_sel))==-counter)) ' / ' num2str(s.basin_sizes(abs(counter)))])
    %     disp(dec2bin(find(ab==counter)-1,n))
    attractor_states_bin=dec2binNum(find(ab==counter)',n);
    basin_states_bin=dec2binNum(find(ab==-counter)',n);
    disp(attractor_states_bin)
        disp(state_names')
        disp('mean:')
        disp(mean(attractor_states_bin,1))
        disp(' ')
    %     disp('mean basin states:')
    %     disp(mean(basin_states_bin))
    %     disp(' ')
    summary=[summary; counter s.attractor_sizes(-counter) s.basin_sizes(-counter) sum(ab(find(ic_sel))==-counter) sum(abs(ab(find(ic_sel)))==-counter)];
    attractor_means=[attractor_means; mean(attractor_states_bin,1)];
end
disp(state_names')

% calculating mean stationary distribution for selected specific ic
attracProb_ic=(summary(:,5)./sum(summary(:,5)))'
pdf_ic=attracProb_ic*attractor_means
pmf_ic=zeros(1,2^n);
for counter=1:length(pmf_ic)
   if (d(counter)==0) %is attractor
      pmf_ic(counter)=attracProb_ic(-ab(counter))/s.attractor_sizes(-ab(counter)); 
   end
end
disp('pmf_ic calculated')