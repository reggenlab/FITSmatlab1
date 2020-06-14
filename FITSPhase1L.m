function FITSPhase1L( varargin ) 

    cargs = varargin ;
    optargin = size(varargin,2);

    passarge ;
    RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
    runID = ceil(100000*rand(1,1) ) + feature('getpid')  
    aRFSR =[60,100];
    maxClusters = 8 ;
    msc = 11 ;
    maxAllowedLevel 
    chunkSize
    
    command = ['awk ''NR == 1'' ', dataName ,' | awk ''{print NF}'' '] ;
    %command
    [status,numberOfSamples] = system(command);
    numberOfSamples = str2num(numberOfSamples);
    %numberOfSamples
    chunks = makeChunks(numberOfSamples, chunkSize) ;  
       
        %tic
        %disp('read')
        %read csv file row consist of sites and column consists of samples
        %dataX = csvread(dataName) ;
        %toc
        %tic
        %disp('norm')
        
    for i = 1:chunks.count
        
        
        %data = dataX(:, chunks.res(i).val);
        ss = size(chunks.res(i).val);
        outfile = strcat(name2save,'_chunk_',num2str(i),'_',num2str(runID),'.txt');
        
        command = ['awk ''{ for(i=',num2str(chunks.res(i).val(ss(1))),'; i<=',num2str(chunks.res(i).val(ss(2))),'; i++) printf "%s",$i OFS; printf ORS}'' ', dataName, ' > ', outfile];
        %command
        system(command);    
        %save(strcat(name2save,'_chunk_',num2str(i),'_',num2str(runID),'.mat'),'data','-v7.3');
    end
    
    %clear dataX;
    disp('chunks creation done');
    for i = 1:chunks.count
        dataX = load(strcat(name2save,'_chunk_',num2str(i),'_',num2str(runID),'.txt'));
        disp(strcat('chunk ', i,' loaded'));
        tic
        %dataX = load(strcat(name2save,'_chunk_',num2str(i),'_',num2str(runID),'.mat'));
        %dataX = dataX.data;
        M = mean(dataX);
        dataX = dataX./(M + 0.00000001);
        data = dataX';
        clear dataX;
        delete(strcat(name2save,'_chunk_',num2str(i),'_',num2str(runID),'.txt'));
        %toc
        disp(strcat('chunk ', i,' normalization done'));
        if fast == 0
            FITSPhase1Start(data,maxClusters,msc,aRFSR,maxAllowedLevel,strcat(name2save,'_fitsL_',num2str(i),'_',num2str(runID)),runID);
        else
            FITSPhase1StartFast(data,maxClusters,msc,aRFSR,maxAllowedLevel,strcat(name2save,'_fitsL_',num2str(i),'_',num2str(runID)),runID, fast);
        end
        indexes = chunks.res(i).val;
        save(strcat(name2save,'_',num2str(i),'_',num2str(runID),'_indexes.mat'),'indexes','-v7.3');
        toc
        disp(strcat('chunk ', i,' completed'));
    end
end

function result = makeChunks(numberOfSamples, chunkSize)
    %numberOfSamples = size(csvread(dataName),2);
    %chunkSize = 1000;
    v = 1:numberOfSamples;
    rv = v;%randsample(v,length(v));
    nonOverlapChunks = ceil(numberOfSamples/chunkSize);
    res = struct;
    start = 1;
    count = 0;
    for i = 1 : nonOverlapChunks
        if i ~= nonOverlapChunks
            res(i).val = rv(start:start+chunkSize-1);
            start = start+chunkSize;
        else
            if i == 1
                res(i).val = rv;
                count = 1;
            else
                if length(rv(start:numberOfSamples)) > 500
                    res(i).val = rv(start:numberOfSamples);
                    count = i;
                else
                    count = i-1;
                    res(i-1).val = rv(start-chunkSize:numberOfSamples);
                end
            end
        end
    end
%     overlapChunks = count;
%     if overlapChunks > 15
%         overlapChunks = 15;
%     end
%     for i =1:overlapChunks
%         res(i+count).val = randsample(v,chunkSize);
%     end
%     count = count + overlapChunks;
    result = struct;
    result.count = count;
    result.res = res;
end
