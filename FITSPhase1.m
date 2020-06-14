function FITSPhase1( varargin ) 

    cargs = varargin ;
    optargin = size(varargin,2);

    passarge ;

    RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
    runID = ceil(100000*rand(1,1) ) + feature('getpid')  
    aRFSR =[60,100];
    maxClusters = 8 ;
    msc = 11 ;
    maxAllowedLevel 
   %tic
    %disp('read')
    %read csv file row consist of sites and column consists of samples
    dataX = load(dataName);%csvread(dataName) ;
    %toc
    %tic
    %disp('norm')
    M = mean(dataX);
    dataX = dataX./(M + 0.00000001);
    dataX = dataX';
    %toc
    if fast == 0
        FITSPhase1Start(dataX,maxClusters,msc,aRFSR,maxAllowedLevel,strcat(name2save,'_complete_',num2str(runID)),runID);
    else
        FITSPhase1StartFast(dataX,maxClusters,msc,aRFSR,maxAllowedLevel,strcat(name2save,'_complete_',num2str(runID)),runID, fast);
    end
end
