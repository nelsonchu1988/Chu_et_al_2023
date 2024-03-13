function Y = modelPredictions(netE,netD,mbq)

Y = [];

% Loop over mini-batches.
while hasdata(mbq)
    X = next(mbq);

    % Forward through encoder.
    Z = predict(netE,X);

    % Forward through dencoder.
    XGenerated = predict(netD,Z);

    % Extract and concatenate predictions.
    Y = cat(4,Y,extractdata(XGenerated));
end

end