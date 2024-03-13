function [Z_cluster, X_cluster, X_cluster_generated, Y_cluster, indexes_cluster] = encoderPredictionsMLP_cluster(netD,netE,mbq,cluster,testset_integers)
    reset(mbq);
    X_cluster = [];
    Z_cluster = [];
    Y_cluster = [];
    indexes_cluster = [];
    X_cluster_generated = [];
    counter = 0;
    % Loop over mini-batches to get the original data and generated data
    while hasdata(mbq)
        counter = counter + 1;
        [X_individual, Y_individual] = next(mbq);
    
        % Forward through encoder.
        %Z_individual = predict(netE,X_individual,Outputs='latentOuput');
        [Z_individual,mu,logSigmaSq] = predict(netE,X_individual);
    
        cluster_index = [];
        for i = 1:size(mu,2)
            mu_individual = mu(:, i);
            mu_x = mu_individual(1);
            mu_y = mu_individual(2);
            cluster_x_start = cluster(1);
            cluster_x_end = cluster(2);
            cluster_y_start = cluster(3);
            cluster_y_end = cluster(4);
            %if mu_individual(1) > cluster(1) && mu_individual(2) > cluster(3) && mu_individual(1) < cluster(2) && mu_individual(2) < cluster(4)
            if mu_x > cluster_x_start && mu_x < cluster_x_end && mu_y > cluster_y_start && mu_y < cluster_y_end
                cluster_index = cat(1, cluster_index, i);
            end
        end
        if isempty(cluster_index)
            continue
        end
        Z_cluster = cat(2,Z_cluster,extractdata(mu(:, cluster_index)));
        X_cluster = cat(2,X_cluster,X_individual(:, cluster_index));
        X_cluster_generated = cat(2,X_cluster_generated,predict(netD, mu(:, cluster_index)));

        Y_individual = extractdata(gather(Y_individual));
        Y_cluster = cat(2,Y_cluster,Y_individual(cluster_index));
        indexes_cluster = cat(2,indexes_cluster, testset_integers((counter - 1)*mbq.MiniBatchSize + cluster_index));
    end
end