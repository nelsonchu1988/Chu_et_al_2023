function [Z, Y] = encoderPredictions(netE,mbq)

Z = [];
Y = [];

% Loop over mini-batches.
while hasdata(mbq)
    [X_individual, Y_individual] = next(mbq);

    % Forward through encoder.
    %Z_individual = predict(netE,X_individual,Outputs='latentOuput');
    [Z_individual,mu,logSigmaSq] = predict(netE,X_individual);
    
    % Extract and concatenate predictions.
    %Z = cat(2,Z,extractdata(Z_individual));
    Z = cat(2,Z,extractdata(mu));

    Y_individual = extractdata(gather(Y_individual));
    
    Y_number = [];
    for col = 1:size(Y_individual, 2)
        Y_number(col) = find(Y_individual(:, col) == 1);
    end
    Y = cat(2,Y,Y_number);

    if size(Y, 2) ~= size(Z, 2)
        szie(Y, 2)
    end
end

end