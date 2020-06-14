dataName = 'xyz' ;
maxAllowedLevel = 4 ;
name2save = 'FITSOutput';
vartrack = zeros(5 ,1) ;
chunkSize = 1000;
fast = 0;
for i = 1:optargin
    checkt = textscan( cargs{i} , '%s' , 'delimiter' , '=') 

    if (strcmp( checkt{1}{1} , 'input' ) == 1)
        dataName = checkt{1}{2} ;
        vartrack(1) = 1 ;
    end

    if (strcmp( checkt{1}{1} , 'output' ) == 1)
       name2save = checkt{1}{2} ;
       vartrack(2) = 1 ;
    end
    
    if (strcmp( checkt{1}{1} , 'maxLevel' ) == 1)
        tempd =  textscan(checkt{1}{2} , '%d')  ;
        maxAllowedLevel  = tempd{1}  ;  
       vartrack(3) = 1 ;
    end

    if (strcmp( checkt{1}{1} , 'chunksize' ) == 1)
       tempd =  textscan(checkt{1}{2} , '%d')  ;
       chunkSize = tempd{1};
       vartrack(4) = 1 ;
    end
    
    
    if (strcmp( checkt{1}{1} , 'fast' ) == 1)
       tempd =  textscan(checkt{1}{2} , '%d')  ;
       fast = tempd{1};
       vartrack(4) = 1 ;
    end
end



disptext = { 
'input=inputCountFile    filename of read-counts csv file, no need to provide genomic location' 
'             '
           
'output=OUPUTfile    filename of imputed matrix to be saved'
'            '

'maxLevel= maxAllowed_level_of_tree  maximum depth of tree  '  
'        '

'chunksize= divide data into small chunks of given size than impute (needed for large datasets default 1000)  '  
'   '

'fast = 0 normal speed(default)  | 2 (fast)  '  
'        '
} ;




if( (vartrack(1) == 0 ) | (vartrack(2) == 0 ) )
    fprintf('%s\n', disptext{:}) ;
end

if(vartrack(1) == 0 )
   disp('errors : No sample data file given' );
   exit(1) ; 
end


if (fast~=0 && fast ~=2)
	disp('fast can take value either 0 or 2');
	exit(1);
end
if(vartrack(2) == 0 )
   disp('warning : No output file name given therefore file save with default name FITSOutput_<??>');
    
end

