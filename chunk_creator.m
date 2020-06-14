function chunk_creator(chunkSize, dataName, name2save)
command = ['awk ''NR == 1'' ', dataName ,' | awk ''{print NF}'' '];
[status,numberOfSamples] = system(command);
numberOfSamples = str2num(numberOfSamples);
chunks = makeChunks(numberOfSamples, chunkSize) ;


for i = 1:chunks.count
	ss = size(chunks.res(i).val);
	outfile = strcat(name2save,'_chunk_',num2str(i),'.txt');

	command = ['awk ''{ for(i=',num2str(chunks.res(i).val(ss(1))),'; i<=',num2str(chunks.res(i).val(ss(2))),'; i++) printf "%s",$i OFS; printf ORS}'' ', dataName, ' > ', outfile];
	system(command);
							         %save(strcat(name2save,'_chunk_',num2str(i),'_',num2str(runID),'.mat'),'data','-v7.3');
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
			if length(rv(start:numberOfSamples)) > (chunkSize/2)
				res(i).val = rv(start:numberOfSamples);
				count = i;
			else
				count = i-1;
				res(i-1).val = rv(start-chunkSize:numberOfSamples);
			end
		end
	end
end
result = struct;
result.count = count;
result.res = res;
end



