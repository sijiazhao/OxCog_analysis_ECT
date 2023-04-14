function myviolin(data, pos, colour)

% args = obj.checkInputs(data, pos, varargin{:});

args.Bandwidth = [];
args.Width = 0.3;

data = data(not(isnan(data)));

hold('on');

% calculate kernel density estimation for the violin
if isempty(data)
    return
end
[density, value] = ksdensity(data, 'bandwidth', args.Bandwidth);
density = density(value >= min(data) & value <= max(data));
value = value(value >= min(data) & value <= max(data));
value(1) = min(data);
value(end) = max(data);

% all data is identical
if min(data) == max(data)
    density = 1;
end

width = args.Width/max(density);

% plot the violin
% patch([pos+density*width pos*ones(size(density))], [value value(end:-1:1)], colour, 'FaceAlpha', 0.3);
patch([pos+density*width pos*ones(size(density))], [value value(end:-1:1)], colour, 'FaceAlpha', 0.3,'EdgeColor','none');
end
