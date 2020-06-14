function FITSPhase1StartFast(unimpData,maxClusters,msc,aRFSR,maxAllowedLevel,name2save,num , fast )
    unimpData=log(unimpData+1.01);
    [ux,uy] = size(unimpData)
    rng shuffle;
    %trees=struct;
    rank = randi([3,7]);
    %tic
    %disp('initial');
    

   K = 30 ;
   if( uy < K  )
       K = uy ;
   end
   [impData,~,~] = svdsecon( unimpData , K);

    if fast ~= 1
       fast = 0;
    end
    %clear dataX;
    RFSR=aRFSR;
    level = 1;
    flag=1;
    clusters = struct;
    records = struct;
    cluster = maxClusters; % help in keeping track of maxclusters allowed in each level
    %class(maxAllowedLevel)
    allowedLevel = 2 ; %randi([2,maxAllowedLevel]); 
    while (level<=allowedLevel) && (flag~=0)
        %tic
        disp('level')
        level
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
                %size(records(level).c(i).val)
                records(level).unimpdata(i).val = unimpData(records(level).c(i).val,:);
                records(level).vec(i).val = find(~any(records(level).unimpdata(i).val, 1));%vector for every data to know column/sites indices having all zero
                records(level).data(i).val = records(level).unimpdata(i).val; %making duplicate matrix
                records(level).data(i).val(:,records(level).vec(i).val) = []; %delete columns having all zero
            
                records(level).unimpdata(i).filename = strcat(name2save,'_unimpdata_level_',num2str(level),'_i_',num2str(i)); 
                records(level).unimpdata(i).val = dumpData(records(level).unimpdata(i).filename,records(level).unimpdata(i).val,1);
                
                [records(level).data(i).val records(level).vec(i).val] = operationOnMatrix(records(level).data(i).val,rank,impData,records(level).c(i).val,records(level).vec(i).val,records(level).vec(i).val, fast);% impute new clustered data
                
                records(level).data(i).filename = strcat(name2save,'_impdata_level_',num2str(level),'_i_',num2str(i)); 
                records(level).data(i).val = dumpData(records(level).data(i).filename,records(level).data(i).val,1);
            end
            clear impData;
            clear unimpData;
            records(level).predictedLabel = 0;
            disp('level_complete');

        else
            try
                disp('level_internal');
        
               records(level).predictedLabel=struct;
               records(level).predictedLabel=clusterPrediction(1,level,clusters,records(level-1).data,RFSR);
            
            
               records(level).c=clusterMaking(1,level,clusters,records(level).predictedLabel);

               records(level).unimpdata = rawDataClusterMaking(1,level,clusters,records(level-1).unimpdata,records(level).c);
               records(level).vec = vecFormation(1,level,clusters,records(level).unimpdata);
               records(level).data = dataFormation(1,level,clusters,records(level).unimpdata,records(level).vec, records(level-1).data);
               [records(level).data records(level).vec] = imputationOperation(rank,1,level,clusters,records(level).data,records(level-1).data,records(level).vec,records(level-1).vec,records(level).c);
                disp('level_internal_complete');
        
            catch
               flag = 0
            end
            records(level).predictedLabel = 0; 
        end
        cluster = ceil(clusters(level).val/2);
        if cluster<2
            cluster=2;
        end
        if (flag~=0)
          level = level + 1;
        end
        %toc
    end
    %input('donot press anything')
    %input('last time not press')
    level=level-1;
    last = level;
    prev = last;
    disp('backtrace')
    while(level>=0)
        %tic
        disp('level')
        level
        if last==level
            records(level).newdata = newdataInitialization(1,level,clusters,records(level).data,uy);
            records(level).vec_negation = vec_negationInitialization(1,level,clusters,records(level).vec,uy);
            records(level).newdata = newdataMaking(1,level,clusters,records(level).newdata,records(level).vec_negation,records(level).data);
        elseif level==0
            pre_final_imputed_new = zeros(ux,uy);
            for i = 1:clusters(1).val
            	pre_final_imputed_new(records(prev).c(i).val,:) = loadData(records(prev).newdata(i).filename, 0, 1);
            end
        else
            records(level).newdata = newdataMerging(1,prev,clusters,records(prev).newdata,records(level).data,records(prev).c,uy);
        end
        prev = level ; 
        level = level - 1;
        %toc
    end
    final_imputed = zeros(ux,uy);
    final_imputed = pre_final_imputed_new;
    kk= dropAllInternalFile(name2save);
    name2save = strcat(name2save,'_',num2str(num),'.mat');
    
    save(name2save,'final_imputed','-v7.3');
    
%toc
end

function  z = istCalc(matrix,rank)
        N1 = numel(matrix);
        %t1 = 1:N1;
        IDX1 =find(matrix>0);
        M1 = opRestriction(N1,IDX1);
        clear N1;
        clear IDX1;
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

function [toreturn updatedvec] = operationOnMatrix(matrix,rank,impmatrix,rows,newVec,oldVec, fast)
    if fast == 1
        toreturn = matrix(rows,:);
        updatedvec = oldVec;
    elseif size(matrix,1)>=rank
        toreturn = istCalc(matrix,rank);
        updatedvec = newVec;
    else
        toreturn = matrix(rows,:);
        updatedvec = oldVec;
    end
end

function predictedLabel = clusterPrediction(startLevel,endLevel,clusters,data,RFSR)
    predictedLabel = struct;
    for i = 1:clusters(startLevel).val
        if startLevel == endLevel-1
        	loaded_data = loadData(data(i).filename, 0, 0);
            predictedLabel(i).val= kclusters(loaded_data,clusters(endLevel).val,RFSR);
            
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
            unimpdata(i).val = loadData(unimpdata(i).filename, 0, 0);
            for m = 1:clusters(endLevel).val
                unimpdata(i).unimpdata(m).val= unimpdata(i).val(c(i).c(m).val,:);
                unimpdata(i).unimpdata(m).filename = strcat(unimpdata(i).filename,'_level_',num2str(endLevel),'_i_',num2str(m)); 
                unimpdata(i).unimpdata(m).val = dumpData(unimpdata(i).unimpdata(m).filename,unimpdata(i).unimpdata(m).val,1);
            end
            unimpdata(i).val = 0;
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
            	loaded_data = loadData(unimpdata(i).unimpdata(m).filename, 0, 0);
                vec(i).vec(m).val= find(~any(loaded_data, 1));
            end
        else
            vec(i).vec=struct;
            vec(i).vec=vecFormation(startLevel+1,endLevel,clusters,unimpdata(i).unimpdata);
        end 
    end
end

%making duplicate matrix
%delete columns having all zero
function data = dataFormation(startLevel,endLevel,clusters,unimpdata,vec,data)
    %data = struct;
    for i = 1:clusters(startLevel).val
        if startLevel == endLevel-1
            for m = 1:clusters(endLevel).val
            	loaded_data = loadData(unimpdata(i).unimpdata(m).filename, 0, 0);
                data(i).data(m).val= loaded_data;
                data(i).data(m).val(:,vec(i).vec(m).val)=[];
                data(i).data(m).filename = strcat(data(i).filename,'_level_',num2str(endLevel),'_i_',num2str(m)); 
                data(i).data(m).val = dumpData(data(i).data(m).filename,data(i).data(m).val,1);
            end
        else
            %data(i).data=struct;
            data(i).data=dataFormation(startLevel+1,endLevel,clusters,unimpdata(i).unimpdata,vec(i).vec, data(i).data);
        end 
    end
end

function [data vec] = imputationOperation(rank,startLevel,endLevel,clusters,data,impdata,vec,oldvec,c)
    for i = 1:clusters(startLevel).val
        if startLevel == endLevel-1
            for m = 1:clusters(endLevel).val
            	loaded_data = loadData(data(i).data(m).filename, 0, 0);
            	loaded_data_imp = loadData(impdata(i).filename, 0, 0);
                fast = 0;
                [data(i).data(m).val vec(i).vec(m).val]= operationOnMatrix(loaded_data,rank,loaded_data_imp,c(i).c(m).val,vec(i).vec(m).val,oldvec(i).val, fast);
                data(i).data(m).filename = strcat(data(i).filename,'_level_',num2str(endLevel),'_i_',num2str(m)); 
                data(i).data(m).val = dumpData(data(i).data(m).filename,data(i).data(m).val,1);
            end
        else
            [data(i).data vec(i).vec]=imputationOperation(rank,startLevel+1,endLevel,clusters,data(i).data,impdata(i).data,vec(i).vec,oldvec(i).vec,c(i).c);
        end 
    end
end

function newdata = newdataInitialization(startLevel,endLevel,clusters,data,dataX_column_size)
    newdata = struct;
    for i = 1:clusters(startLevel).val
        if startLevel == endLevel-1
            for m = 1:clusters(endLevel).val
            	row_size = loadData(data(i).data(m).filename, 1, 0);
                newdata(i).newdata(m).val= zeros(row_size,dataX_column_size);
                newdata(i).newdata(m).filename = strcat(data(i).data(m).filename,'_newdata'); 
                newdata(i).newdata(m).val = dumpData(newdata(i).newdata(m).filename,newdata(i).newdata(m).val,1);
            end
        else
            newdata(i).newdata=struct;
            newdata(i).newdata=newdataInitialization(startLevel+1,endLevel,clusters,data(i).data,dataX_column_size);
        end 
    end
    
end 

function vec_negation = vec_negationInitialization(startLevel,endLevel,clusters,vec,dataX_column_size)
    vec_negation = struct;
    for i = 1:clusters(startLevel).val
        if startLevel == endLevel-1
            for m = 1:clusters(endLevel).val
                vec_negation(i).vec_negation(m).val=  1:dataX_column_size;
                vec_negation(i).vec_negation(m).val(vec(i).vec(m).val)=[];
                     
            end
        else
            vec_negation(i).vec_negation=struct;
            vec_negation(i).vec_negation=vec_negationInitialization(startLevel+1,endLevel,clusters,vec(i).vec,dataX_column_size);
        end 
    end
    
end 

function newdata = newdataMaking(startLevel,endLevel,clusters,newdata,vec_negation,data)
    for i = 1:clusters(startLevel).val
        if startLevel == endLevel-1
            for m = 1:clusters(endLevel).val
            	loaded_data = loadData(data(i).data(m).filename, 0, 0);
            	newdata(i).newdata(m).val = loadData(newdata(i).newdata(m).filename, 0, 0);
                newdata(i).newdata(m).val(:,vec_negation(i).vec_negation(m).val) = loaded_data;
                newdata(i).newdata(m).val = dumpData(newdata(i).newdata(m).filename,newdata(i).newdata(m).val,1);
            end
        else
            newdata(i).newdata=newdataMaking(startLevel+1,endLevel,clusters,newdata(i).newdata,vec_negation(i).vec_negation,data(i).data);
        end 
    end
end

function newdata = newdataMerging(startLevel,endLevel,clusters,newdata,data,c,dataX_column_size)
    for i = 1:clusters(startLevel).val
        if startLevel == endLevel-1
        	row_size = loadData(data(i).filename, 1, 0);
            newdata(i).val = zeros(row_size,dataX_column_size);
            for m = 1:clusters(endLevel).val
            	loaded_data = loadData(newdata(i).newdata(m).filename, 0, 1); %matrix deleted here
                newdata(i).val(transpose(c(i).c(m).val),:) = loaded_data;
            end
            newdata(i).filename = strcat(data(i).filename,'_newdata'); 
            newdata(i).val = dumpData(newdata(i).filename,newdata(i).val,1);
        else
            newdata(i).newdata=newdataMerging(startLevel+1,endLevel,clusters,newdata(i).newdata,data(i).data,c(i).c,dataX_column_size);
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

function y = dumpData(name2save, internalData, isReturnZero)
    save(strcat(name2save,'_internal_chdsffuegfcalwc12123v2rjvgv32hyvfh32yrh.mat'),'internalData','-v7.3');
    internal_data_size = size(internalData);
    save(strcat(name2save,'_internal_size_chdsffuegfcalwc12123v2rjvgv32hyvfh32yrh.mat'),'internal_data_size','-v7.3');
    if isReturnZero == 1
        y = 0;
    else
        y = internalData;
    end
end

function y = loadData(name2load, isReturnSize, isDeleteFile)
    if isReturnSize ~= 0
        internal_obj = load(strcat(name2load,'_internal_size_chdsffuegfcalwc12123v2rjvgv32hyvfh32yrh.mat'));
        y = internal_obj.internal_data_size(isReturnSize);
    else
        internal_obj = load(strcat(name2load,'_internal_chdsffuegfcalwc12123v2rjvgv32hyvfh32yrh.mat'));
        y = internal_obj.internalData;
    end
    
    if  isDeleteFile == 1
        %disp('deleted file')
        %strcat(name2load,'_internal_chdsffuegfcalwc12123v2rjvgv32hyvfh32yrh.mat')
        delete(strcat(name2load,'_internal_chdsffuegfcalwc12123v2rjvgv32hyvfh32yrh.mat'));
        delete(strcat(name2load,'_internal_size_chdsffuegfcalwc12123v2rjvgv32hyvfh32yrh.mat'));
    end
end

function dropDone = dropAllInternalFile(name2save)
    dropDone = 0;
    files = dir(strcat(name2save,'*internal_chdsffuegfcalwc12123v2rjvgv32hyvfh32yrh.mat'));
    for t = 1:size(files,1)
        delete(strcat(files(t).folder,'/',files(t).name));
    end

    files = dir(strcat(name2save,'*internal_size_chdsffuegfcalwc12123v2rjvgv32hyvfh32yrh.mat'));
    for t = 1:size(files,1)
        delete(strcat(files(t).folder,'/',files(t).name));
    end

end
