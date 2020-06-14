"# FITS Matlab version" 

You have to download matlab code in your local machine/server. For execution you have to pass filename of read-counts csv as a input file. Read-count csv file does not contain header and genomic location i.e. it consist of only data on which imputation is going to perform. Row represents sites and column reprent samples/cells in csv file.

"# For small dataset"

```matlab
FITSPhase1 input='<csv file name>'
e.g.
FITSPhase1 input='sce5_raw.csv'
```
'sce5_raw.csv' consist of epignome data corresponding to five cell type.

Other optional input parameter you can pass in phase1 

```matlab
FITSPhase1 input='<csv file name>' output='<name to save imputed file>' maxLevel =<Depth upto which tree will grow> fast=<0/2>
```

By default maxLevel set to 4, fast =0  and output set to 'FITSOutput'.
You can run FITSPhase1 parallely in background using.
Parameter fast=2 run FITS faster as compared to fast=0 but it limits you from generic to specific default parameters. 
Parameter output='folderpath/nametosave' name to save do not end with any extension e.g. if you have to save file in abc folder with name sce5.csv then just pass output='abc/sce5' 

```bash
nohup matlab -nodisplay -nosplash -r "try FITSPhase1 input='<csv file name>'; catch; end; quit" > <name>.txt &
```
You can create n number of imputed matrix generated through phase1. Each run will generate imputed matrix.

Once Phase1 is over then run Phase2 to generate final imputed matrix based on matrix received as output from Phase1.

```matlab
FITSPhase2 input='<csv file name>'
e.g.
FITSPhase2 input='sce5_raw.csv'
```
You have to pass same input file as you passed in Phase1. Don't worry Phase2 takes only one minute to generate final output :)

Other optional input parameter you can pass in phase2 

```matlab
FITSPhase2 input='<csv file name>' output='<name to save imputed file same as Phase1>' k =<topk correlated matrix feature/sample value use for final imputation> feature=<1/0 takes values either 1 or 0>  remove=<1/0> saveformat=<'csv'/'tab' > 
```
Default value for feature is zero. At value 0 phase2 will compute correlation among samples/cell (preffered) while value 1 will compute correlation features/sites wise.
saveformat takes 'csv'/'tab' parameters to save final imputed matrix in format you want, default (csv). If 'tab' pass then you get tab separated ouput file saved with extension .txt
remove takes values either o(default) or 1 (if want to delete all intermediate files generated (recommended)) 



"# For large dataset"

```matlab
FITSPhase1L input='<csv file name>'
e.g.
FITSPhase1L input='sce5_raw.csv'
```
'sce5_raw.csv' consist of epignome data corresponding to five cell type.

Other optional input parameter you can pass in phase1L 

```matlab
FITSPhase1L input='<csv file name>' output='<name to save imputed file>' maxLevel =<Depth upto which tree will grow> fast=<0/2> chunksize=<default 1000, size in which you want to divide samples eg. 1500>
```
more about chunksize: 
chunksize will make your big data file into small. If your file consist of 10,200 samples and you passed chunk size value =1000 then it will create 10 files in total of 1000 non overlapping samples each with exception (1 file will consist of 1200 samples). But if there are 10,900 samples instead of 10,200 then 11 files will be created with exception (1 file consist of only 900 samples). 

You can run n number of FITSPhase1L process parallely depending upon your system computation power
> keep in mind, to pass same output name in all parallel process for same data

Once Phase1L is over you can run Phase2L

```matlab
FITSPhase2L input='<csv file name>'
e.g.
FITSPhase2L input='sce5_raw.csv'
```
You have to pass same input file as you passed in Phase1. Don't worry Phase2 takes only one minute to generate final output :)

Other optional input parameter you can pass in phase2 

```matlab
FITSPhase2L input='<csv file name>' output='<name to save imputed file same as Phase1L>' k =<topk correlated matrix feature/sample value use for final imputation>  remove=<1/0> saveformat=<'csv'/'tab'> 
```

**How To Decide, When to use FITSPhase1 or FITSPhase1L?**
It completely depends on your system computation power. We recommend to go for FITSPhase1 if you have computation power as it does not break data into chunks. hence will impute on complete data whereas FITSPhase1L will break data into chunks and then impute (although you have not to worry how chunks are created as in end you get merged data). 

**How To Decide, When to use FITSPhase2 or FITSPhase2L?**
Simple, if you computed imputation through FITSPhase1 then use FITSPhase2 otherwise FITSPhase2L
