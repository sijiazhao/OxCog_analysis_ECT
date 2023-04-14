function y = bootstrap2(data,n)
[NS, t] = size(data); % NS = number of sample, t = time sample

B = 5000; %number of iterations
% p = repmat(1:NS,B,1);
% p = p(reshape(randperm(B*NS),B,NS));

y = zeros(t,B);

for b=1:B
    p = datasample(1:NS,n,'Replace',true);
%     length(unique(p))
    y(:,b) = nanmean(data(p,:));
end

end
