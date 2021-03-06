function c = PTPalette(n)
% Colors from Paul Tol's qualitative color palette
% http://www.sron.nl/~pault/colourschemes.pdf

if (n > 12)
    disp('Only 12 colors defined in palette');
    c = [];
    return;
end

colors = [...
    51 34 136; ...
    204 102 119; ...
    221 204 119; ...
    17 119 51; ...
    136 204 238; ...
    170 68 153; ...
    68 170 153; ...
    153 153 51; ...
    136 34 85; ...
    102 17 0; ...
    102 153 204; ...
    170 68 102; ...
    ]./255;

c = colors(1:n,:);
if (n <=4)
    c(1,:) = [68 119 170]./255;
end


end

