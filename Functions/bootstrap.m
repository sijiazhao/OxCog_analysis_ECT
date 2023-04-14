function y = bootstrap(data)
[NS, t] = size(data); % NS = number of sample, t = time sample

meandata = mean(data);

B = 10000; %number of iterations
p = repmat(1:NS,B,1);
p = p(reshape(randperm(B*NS),B,NS));

y = zeros(length(data),B);

for b=1:B
    y(:,b) = mean(data(p(b,:),:));
end

end
