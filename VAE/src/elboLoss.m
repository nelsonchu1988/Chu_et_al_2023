function loss = elboLoss(Y,T,mu,logSigmaSq)

% Reconstruction loss.
%reconstructionLoss = crossentropy(Y,T,'TargetCategories','independent');
reconstructionLoss = mse(Y,T);

% KL divergence.
KL = -0.5 * sum(1 + logSigmaSq - mu.^2 - exp(logSigmaSq),1);
KL = mean(KL);

% Combined loss.
loss = reconstructionLoss + KL;

%imshow(extractdata(squeeze(Y(:,:,1,1))))
%imshow(extractdata(squeeze(T(:,:,1,1))))
end