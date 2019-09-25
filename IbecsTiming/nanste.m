function [ e ] = nanste(M)
%
s = size(M,1)
e = nanstd(M)./sqrt(s-1)

end

