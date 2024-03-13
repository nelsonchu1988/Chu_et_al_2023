function [X, Y] = preprocessMiniBatch(XCell, YCell)

% Concatenate.
X = cat(4,XCell{:});

% Extract label data from cell and concatenate.
Y = cat(2,YCell{:});

% One-hot encode labels.
%Y = onehotencode(Y,1);

end