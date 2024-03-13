function [X, Y] = load_lightSheet_data()
    if isfile([pwd filesep 'data' filesep 'lightSheetData.mat'])
        data = load([pwd filesep 'data' filesep 'lightSheetData.mat']);
        X = data.X;
        Y = data.Y;
        return
    end
    r = 0:10:200;
    dim = 'CXCL12';
    subs = {'OP1', 'OP2','OP3', 'YP3', 'YP4', 'YP5'};
    X = [];
    Y = "";
    counter = 1;
    for i = 1:length(subs)
        for j = 1:length(r)
            resultFile = [pwd filesep 'data' filesep subs{i} '_dist_' dim '_' num2str(r(j)) '.mat'];
            data = load(resultFile);
            
            % construct input data X
            X(1, j, 1, counter:counter + length(data.distances) - 1) = data.distances.cxcl12Count;
            X(2, j, 1, counter:counter + length(data.distances) - 1) = data.distances.boneCount;
            X(3, j, 1, counter:counter + length(data.distances) - 1) = data.distances.cd31Count;
            X(4, j, 1, counter:counter + length(data.distances) - 1) = data.distances.yoproCount;
            X(5, j, 1, counter:counter + length(data.distances) - 1) = data.distances.periCount;
            
            % fill zero 
            X(6:28, j, 1, counter:counter + length(data.distances) - 1) = 0;

            % construct labels Y
            sub = subs{i};
            Y(1, counter:counter + length(data.distances) - 1) = sub(1);
            Y(2, counter:counter + length(data.distances) - 1) = sub;
            Y = categorical(Y);
        end
        counter = counter + length(data.distances);
        X(:, 22:28, 1, counter:counter + length(data.distances) - 1) = 0;
    end
    save([pwd filesep 'data' filesep 'lightSheetData.mat'], "X", "Y", "-v7.3");
end