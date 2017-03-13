function [outv] = downsamplebin(invec,amount)
%discards last few time points that don't fit in a bin anymore


nbins=floor(length(invec)/amount);
outv = sum(reshape(invec(1:nbins*amount),amount,nbins) ,1);
