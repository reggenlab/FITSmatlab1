function FITSPhase2(varargin)
    
    cargs = varargin ;
    optargin = size(varargin,2);

    passarge2 ;
    files = dir(strcat(name2save,'_complete_*.mat'));
    dataX = load(dataName) ;
    M = mean(dataX);
    dataX = dataX./(M + 0.00000001);
    dataX = dataX';
    dataX=log(dataX+1.01);
    trees=struct;
    start=1;
    for t = 1:size(files,1)
      if not(contains(files(t).name , '_chdsffuegfcalwc12123v2rjvgv32hyvfh32yrh.mat'))
        try
           obj = load(strcat(files(t).folder ,'/', files(t).name));
           trees(start).val = obj.final_imputed;
           %size(trees(start).val);
           start = start+1;
        catch
           %numOfTrees = numOfTrees-1;
           disp(strcat(files(t).folder ,'/', files(t).name));
           disp('This is either corrupted or not exist');
        end
      end
    end
    [row, col] = size(dataX);
    if (topk > start-1)
       topk = start-1;
    end
    maxCorrelated(dataX,trees,start-1,topk,colWise,name2save,csv_or_tab);
    %save(strcat(name2save,'_.mat'),'final_imputed','-v7.3');
    if isDelete == 1
        uu = dropAllInternalFile(name2save)
    end
end

function res = maxCorrelated(mOriginal,mtree,count,topk,colWise,name2save,csv_or_tab)
    if colWise==1
        final_imputed = maxCorrelatedCol(mOriginal,mtree,count,topk);
        %save(strcat(name2save,'.mat'),'final_imputed','-v7.3');
        
    else
        final_imputed = maxCorrelatedRow(mOriginal,mtree,count,topk);
        %save(strcat(name2save,'.mat'),'final_imputed','-v7.3');
    end

    if csv_or_tab ~= 'tab'
    	csvwrite(strcat(name2save,'.csv'),final_imputed');
    else
    	writematrix(final_imputed',strcat(name2save,'.txt'),'Delimiter','tab');
    end
end

% correlation between features
function res = maxCorrelatedCol(mOriginal,mtree,count,topk)
  [row, col] = size(mOriginal);
    res = zeros(row,col);
    for i = 1 : col
        corrAll = [];
        for j = 1 : count
           cor = corr(mtree(j).val(:,i),mOriginal(:,i),'Type','Spearman');
           corrAll = [corrAll; cor];
        end
        [~, c_order] = sort(corrAll,'descend');
        newMatrix = zeros(row,topk);
        for j = 1:topk
            newMatrix(:,j) = mtree(c_order(j)).val(:,i);
        end
        res(:,i) = max(newMatrix')';
    end
end

function dropDone = dropAllInternalFile(name2save)
    dropDone = 0;
    files = dir(strcat(name2save,'_complete_*.mat'));
    for t = 1:size(files,1)
        disp('deleting')
        strcat(files(t).folder,'/',files(t).name)
        delete(strcat(files(t).folder,'/',files(t).name));
    end
end
% %correlation between samples
function res = maxCorrelatedRow(mOriginal,mtree,count,topk)
    [row, col] = size(mOriginal);
    res = zeros(row,col);
    for i = 1 : row
        corrAll = [];
        for j = 1 : count
           cor = corr(mtree(j).val(i,:)',mOriginal(i,:)','Type','Spearman');
           corrAll = [corrAll; cor];
        end
        [~, c_order] = sort(corrAll,'descend');
        newMatrix = zeros(topk,col);
        for j = 1:topk
            newMatrix(j,:) = mtree(c_order(j)).val(i,:);
        end
        res(i,:) = max(newMatrix);
    end
end
