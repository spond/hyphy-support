LIKELIHOOD_FUNCTION_OUTPUT = 1; 
DataSet raw_data = ReadDataFile("data/hiv.seq");
DataSetFilter filt_data = CreateFilter(raw_data,1);
HarvestFrequencies(F4,filt_data,1,1,1); // F4 contains empirical nucleotide frequencies 

F1 = {{1.}, {1.}, {1.}, {1.}} * (0.25); // F1 is simply a vector of ones.

global freqA:=F4__[0];
global freqC:=F4__[1];
global freqG:=F4__[2];
global freqT:=F4__[3];
global k=1;
global trvs=1;

// Regular matrix
HKY85 = {{,trvs,k*trvs,trvs}
{trvs,,trvs,k*trvs}
{k*trvs,trvs,,trvs}
{trvs,k*trvs,trvs,}};

// Matrix which ALREADY contains equilibrium frequencies
HKY85F = {{,trvs*freqC,k*trvs*freqG,trvs*freqT}
{trvs*freqA,,trvs*freqG,k*trvs*freqT}
{k*trvs*freqA,trvs*freqC,,trvs*freqT}
{trvs*freqA,k*trvs*freqC,trvs*freqG,}};

LIKELIHOOD_FUNCTION_OUTPUT = 2;

Model MyModel = (HKY85, F4, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree Tree01 = DATAFILE_TREE;
LikelihoodFunction LikFn = (filt_data, Tree01);
Optimize (paramValues, LikFn);
fprintf (stdout, LikFn);

Model MyModel = (HKY85F, F4, 0);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree Tree01 = DATAFILE_TREE;
LikelihoodFunction LikFn = (filt_data, Tree01);
Optimize (paramValues, LikFn);
fprintf (stdout, LikFn);

Model MyModel = (HKY85F, F1, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree Tree01 = DATAFILE_TREE;
LikelihoodFunction LikFn = (filt_data, Tree01);
Optimize (paramValues, LikFn);
fprintf (stdout, LikFn);

Model MyModel = (HKY85F, F1, 0);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree Tree01 = DATAFILE_TREE;
LikelihoodFunction LikFn = (filt_data, Tree01);
Optimize (paramValues, LikFn);
fprintf (stdout, LikFn);