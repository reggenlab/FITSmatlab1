function FITSPhase1Start(unimpData,maxClusters,msc,aRFSR,maxAllowedLevel,name2save,num)
   
    %obj = load(dataName);
    %dataX = obj.dataX;
    %size(dataX)
    unimpData=log(unimpData+1.01);
    %unimpData = dataX;
    %clear dataX;
    rng shuffle;
    %trees=struct;
        rank = randi([3,6]);
        %tic
        %disp('initial')
        Xrec = istCalc(unimpData,rank);
        final_imputed = Xrec;
        %toc
        save(strcat(name2save,'_r_',num2str(rank),'.mat'),'final_imputed','-v7.3');
        impData = Xrec;
        clear final_imputed;
        clear Xrec;
        %clear dataX;
        RFSR=aRFSR;
        level = 1;
        flag=1;
        clusters = struct;
        records = struct;
        cluster = maxClusters; % help in keeping track of maxclusters allowed in each level
        %class(maxAllowedLevel)
        allowedLevel = randi([2,maxAllowedLevel]); 
        while (level<=allowedLevel) && (flag~=0)
            %tic
            %disp('level')
            %level
            clusters(level).val =  clusterCount(cluster);
            if level==1
                while 1
                    records(level).predictedLabel=kclusters(impData,clusters(level).val,RFSR);
                    %%Data Check function should be implemented here
                    flag = checkClustersMinCount(records(level).predictedLabel,clusters(level).val,msc);
                    if flag~=0
                        break;
                    else
                        clusters(level).val=clusters(level).val-1;
                    end
                end
                for i = 1:clusters(level).val
                    records(level).c(i).val = find(records(level).predictedLabel==i);
                    records(level).unimpdata(i).val = unimpData(records(level).c(i).val,:);
                end


                %vector for every data to know column/sites indices having all zero
                %making duplicate matrix
                %delete columns having all zero
                for i = 1:clusters(level).val
                    records(level).vec(i).val = find(~any(records(level).unimpdata(i).val, 1));
                    records(level).data(i).val = records(level).unimpdata(i).val;
                    records(level).data(i).val(:,records(level).vec(i).val) = [];
                end

                % impute new clustered data
                for i = 1:clusters(level).val
                    [records(level).data(i).val records(level).vec(i).val] = operationOnMatrix(records(level).data(i).val,rank,impData,records(level).c(i).val,records(level).vec(i).val,records(level).vec(i).val);
                end
            else
                try
                   records(level).predictedLabel=struct;
                   records(level).predictedLabel=clusterPrediction(1,level,clusters,records(level-1).data,RFSR);
                
                
                   records(level).c=clusterMaking(1,level,clusters,records(level).predictedLabel);

                   records(level).unimpdata = rawDataClusterMaking(1,level,clusters,records(level-1).unimpdata,records(level).c);
                   records(level).vec = vecFormation(1,level,clusters,records(level).unimpdata);
                   records(level).data = dataFormation(1,level,clusters,records(level).unimpdata,records(level).vec);
                   [records(level).data records(level).vec] = imputationOperation(rank,1,level,clusters,records(level).data,records(level-1).data,records(level).vec,records(level-1).vec,records(level).c);
                catch
                   flag = 0
                end 
            end
            cluster = ceil(clusters(level).val/2);
            if cluster<2
                cluster=2;
            end
            if (flag~=0)
              level = level + 1
            end
            %toc
        end
        level=level-1;
        last = level;
        prev = last;
        %disp('backtrace')
        while(level>=0)
            %tic
            %disp('level')
            %level
            if last==level
                records(level).newdata = newdataInitialization(1,level,clusters,records(level).data,unimpData);
                records(level).vec_negation = vec_negationInitialization(1,level,clusters,records(level).vec,unimpData);
                records(level).newdata = newdataMaking(1,level,clusters,records(level).newdata,records(level).vec_negation,records(level).data);
            elseif level==0
                pre_final_imputed_new = zeros(size(unimpData));
                for i = 1:clusters(1).val
                    pre_final_imputed_new(records(prev).c(i).val,:) = records(prev).newdata(i).val;
                end
            else
                records(level).newdata = newdataMerging(1,prev,clusters,records(prev).newdata,records(level).data,records(prev).c,unimpData);
            end
            prev = level ; 
            level = level - 1;
            %toc
        end
        final_imputed = zeros(size(unimpData));
        final_imputed = pre_final_imputed_new;
        name2save = strcat(name2save,'_',num2str(num),'.mat');
        save(name2save,'final_imputed','-v6');
%toc
end

function  z = istCalc(matrix,rank)
        N1 = numel(matrix);
        %t1 = 1:N1;
        IDX1 =find(matrix>0);
        M1 = opRestriction(N1,IDX1);
        clear N1
        clear IDX1
        %y1 = M1(matrix(:),1);
        z = IST_eMC( M1(matrix(:),1),M1,size(matrix), rank); 
end
function plabels = kclusters(matrix,k,RFSR)
    random_vec = randperm(size(matrix,2));
    r=randi(RFSR)/100;
    val = size(matrix,2)-int32(size(matrix,2)*r);
    matrix(:,random_vec(1:val)) = [];
    %[matrix, mapping] = compute_mapping(matrix, 'PCA',8)
    [m,n] = size(matrix);
    if ~(10 <= m && 10 <= n)
       disp('');
    else
       [U,S,V] = svdsecon(matrix,10);
       [matrix, mapping] = compute_mapping(U, 'tSNE',randi([3,6]));
    end
    %[U,S,V] = svdsecon(matrix,10);
    %disp('size u')
    %size(U)
    %[matrix, mapping] = compute_mapping(U, 'tSNE',randi([3,6]));
    
    %size(matrix,1)
    loc=randperm(size(matrix,1),k);
    init =matrix(loc,:);
    plabels=kmeans(matrix,k,'MaxIter',1000,'Start',init);
end 
function v = clusterCount(maxk)
    val=randperm(maxk);
    v=val(1);
    if v<2
        v=2;
    end
end

function [toreturn updatedvec] = operationOnMatrix(matrix,rank,impmatrix,rows,newVec,oldVec)
    if size(matrix,1)>=rank
        toreturn = istCalc(matrix,rank);
        updatedvec = newVec;
    else
        toreturn = impmatrix(rows,:);
        updatedvec = oldVec;
    end
end

function predictedLabel = clusterPrediction(startLevel,endLevel,clusters,data,RFSR)
    predictedLabel = struct;
    for i = 1:clusters(startLevel).val
        if startLevel == endLevel-1
            predictedLabel(i).val= kclusters(data(i).val,clusters(endLevel).val,RFSR);
            
        else
            predictedLabel(i).predictedLabel=struct;
            predictedLabel(i).predictedLabel=clusterPrediction(startLevel+1,endLevel,clusters,data(i).data,RFSR);
        end 
    end
end

function c = clusterMaking(startLevel,endLevel,clusters,predictedLabel)
    c = struct;
    for i = 1:clusters(startLevel).val
        if startLevel == endLevel-1
            for m = 1:clusters(endLevel).val
                c(i).c(m).val= find(predictedLabel(i).val==m);
            end
        else
            c(i).c=struct;
            c(i).c=clusterMaking(startLevel+1,endLevel,clusters,predictedLabel(i).predictedLabel);
        end 
    end
end

%original data and their labels in 8 clusters
function unimpdata = rawDataClusterMaking(startLevel,endLevel,clusters,unimpdata,c)
    for i = 1:clusters(startLevel).val
        if startLevel == endLevel-1
            for m = 1:clusters(endLevel).val
                unimpdata(i).unimpdata(m).val= unimpdata(i).val(c(i).c(m).val,:);
            end
        else
            unimpdata(i).unimpdata=rawDataClusterMaking(startLevel+1,endLevel,clusters,unimpdata(i).unimpdata,c(i).c);
        end 
    end
end

%vector for every data to know column/sites indices having all zero
function vec = vecFormation(startLevel,endLevel,clusters,unimpdata)
    vec = struct;
    for i = 1:clusters(startLevel).val
        if startLevel == endLevel-1
            for m = 1:clusters(endLevel).val
                vec(i).vec(m).val= find(~any(unimpdata(i).unimpdata(m).val, 1));
            end
        else
            vec(i).vec=struct;
            vec(i).vec=vecFormation(startLevel+1,endLevel,clusters,unimpdata(i).unimpdata);
        end 
    end
end

%making duplicate matrix
%delete columns having all zero
function data = dataFormation(startLevel,endLevel,clusters,unimpdata,vec)
    data = struct;
    for i = 1:clusters(startLevel).val
        if startLevel == endLevel-1
            for m = 1:clusters(endLevel).val
                data(i).data(m).val= unimpdata(i).unimpdata(m).val;
                data(i).data(m).val(:,vec(i).vec(m).val)=[];
            end
        else
            data(i).data=struct;
            data(i).data=dataFormation(startLevel+1,endLevel,clusters,unimpdata(i).unimpdata,vec(i).vec);
        end 
    end
end

function [data vec] = imputationOperation(rank,startLevel,endLevel,clusters,data,impdata,vec,oldvec,c)
    for i = 1:clusters(startLevel).val
        if startLevel == endLevel-1
            for m = 1:clusters(endLevel).val
                [data(i).data(m).val vec(i).vec(m).val]= operationOnMatrix(data(i).data(m).val,rank,impdata(i).val,c(i).c(m).val,vec(i).vec(m).val,oldvec(i).val);
            end
        else
            [data(i).data vec(i).vec]=imputationOperation(rank,startLevel+1,endLevel,clusters,data(i).data,impdata(i).data,vec(i).vec,oldvec(i).vec,c(i).c);
        end 
    end
end

function newdata = newdataInitialization(startLevel,endLevel,clusters,data,dataX)
    newdata = struct;
    for i = 1:clusters(startLevel).val
        if startLevel == endLevel-1
            for m = 1:clusters(endLevel).val
                newdata(i).newdata(m).val= zeros(size(data(i).data(m).val,1),size(dataX,2));
            end
        else
            newdata(i).newdata=struct;
            newdata(i).newdata=newdataInitialization(startLevel+1,endLevel,clusters,data(i).data,dataX);
        end 
    end
    
end 

function vec_negation = vec_negationInitialization(startLevel,endLevel,clusters,vec,dataX)
    vec_negation = struct;
    for i = 1:clusters(startLevel).val
        if startLevel == endLevel-1
            for m = 1:clusters(endLevel).val
                vec_negation(i).vec_negation(m).val=  1:size(dataX,2);
                vec_negation(i).vec_negation(m).val(vec(i).vec(m).val)=[];
                     
            end
        else
            vec_negation(i).vec_negation=struct;
            vec_negation(i).vec_negation=vec_negationInitialization(startLevel+1,endLevel,clusters,vec(i).vec,dataX);
        end 
    end
    
end 

function newdata = newdataMaking(startLevel,endLevel,clusters,newdata,vec_negation,data)
    for i = 1:clusters(startLevel).val
        if startLevel == endLevel-1
            for m = 1:clusters(endLevel).val
                newdata(i).newdata(m).val(:,vec_negation(i).vec_negation(m).val) =data(i).data(m).val;
            end
        else
            newdata(i).newdata=newdataMaking(startLevel+1,endLevel,clusters,newdata(i).newdata,vec_negation(i).vec_negation,data(i).data);
        end 
    end
end

function newdata = newdataMerging(startLevel,endLevel,clusters,newdata,data,c,dataX)
    for i = 1:clusters(startLevel).val
        if startLevel == endLevel-1
            newdata(i).val = zeros(size(data(i).val,1),size(dataX,2));
            for m = 1:clusters(endLevel).val
                newdata(i).val(transpose(c(i).c(m).val),:) = newdata(i).newdata(m).val;
            end
        else
            newdata(i).newdata=newdataMerging(startLevel+1,endLevel,clusters,newdata(i).newdata,data(i).data,c(i).c,dataX);
        end 
    end
end

function res = maxCorrelated(mOriginal,mtree,count)
    [row, col] = size(mOriginal);
    res = zeros(row,col);
    for i = 1 : col
        pos = -1;
        maxcorr = -1;
        for j = 1 : count
            if pos==-1
                pos = j;
                maxcorr = corr(mtree(j).val(:,i),mOriginal(:,i));
            else
                cor = corr(mtree(j).val(:,i),mOriginal(:,i));
                if maxcorr<cor
                    pos = j;
                    maxcorr=cor;
                end
            end
        end
        res(:,i) = mtree(pos).val(:,i);
    end
end
function flag  = checkClustersMinCount(predicted_labels,clusters,msc)
    flag = 1 ;          
    for i = 1:clusters
        val = find(predicted_labels==i);
        if size(val)<msc
            flag = 0;
        end
    end

end
