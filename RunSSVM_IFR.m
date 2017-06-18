This example is run on the Duke University influenza data set


[ Data.TrainData Data.TrainLabels Data.TestData Data.TestLabels] = GetDukeInfluenzaData( 1 );

Data.TOL = 0.001; %interior point tolerance
Data.method = 'ReducedKKTy'
Data.C = 1;%should optimize this SVM paramter in general

NumIterations = 2;%number of time to recursively remove features.
for kk = 1:NumIterations%iterative removal loop
    kk%run one classification
    Data.Model{kk} = SSVM_Train(Data.TrainData, Data.TrainLabels, Data.C, Data.TOL, Data.method);
    %determine accuracy
    Data.ACC{kk} = get_accuracy(Data.TrainData, Data.TrainLabels, Data.Model{kk});
    Data.ACC{kk}.acc
    Data.ACCtest{kk} = get_accuracy(Data.TestData, Data.TestLabels, Data.Model{kk});
    Data.ACCtest{kk}.acc
    [Data.aa{kk} Data.bb{kk}] = sort( abs(Data.Model{kk}.wgt) , 'descend');
    count = 1
%     while Data.aa{kk}(count) >  0.000001
%         count = count + 1;
%     end
    while Data.aa{kk}(count)/Data.aa{kk}(count+1) <  100%secret sauce
        count = count + 1;
    end
    
    count;
    Data.numkept{kk} = count - 1;%how many features were kept this iteration.
    Data.TrainData(:, Data.bb{kk}(1:count-1))  = 0; %%%careful, work with raw indices
    %Genes{kk} = Fluz.gene1(bb{kk}(1:count-1),2);%get the gene names.
    if Data.ACCtest{kk}.acc < .75
        break
    end
end

figure
hold on
for kk = 1:size(Data.ACCtest{kk},2)
    plot(Data.ACCtest{kk}.acc,'.')
    plot(Data.ACC{kk}.acc,'o')
end

