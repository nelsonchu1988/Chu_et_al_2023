function [Z, Y] = encoderPredictionsMLP(netE,mbq)

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
    
    Y = cat(2,Y,Y_individual);
end

end