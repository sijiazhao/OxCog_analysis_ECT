function z = nanzscore_sz(x)
z = (x - nanmean(x))./nanstd(x);
end
