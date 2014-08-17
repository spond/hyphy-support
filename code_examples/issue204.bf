LoadFunctionLibrary("ProbabilityDistributions");
siteRates = {1, 100}["sampleFromGamma(0.5,0.5)"];

/* set this variable to 1 to get HyPhy to 
automatically set model parameter values (like t below)
to match the branch lengths supplied with the tree */

AUTOMATICALLY_CONVERT_BRANCH_LENGTHS = 1;

/* set siteScale to one at first to get the 
branch lengths right */

global siteScale := 1;


HKYM = {4, 4};
HKYs = {{    *, 0.001,   0.2, 0.001}
        {0.001,     *, 0.001,   0.2}
        {  0.2, 0.001,     *, 0.001}
        {0.001,   0.2, 0.001,     *}};

for(i=0; i<4; i += 1){
  for(j=0; j<4; j += 1){
    if(i != j){
      HKYM[i][j] := HKYs[i__][j__] * siteScale * t;
    }
  }
}

eqf = {{0.4,0.1,0.2,0.3}};
Model HKY = (HKYM, eqf, 1);
Tree T = (Chimpanzee:0.0561983,Human:0.0411321,(Gorilla:0.0630279,(Gibbon:0.132441,Orangutan:0.0979926)Node6:0.053019)Node4:0.0143057);

fprintf (stdout, "Using this tree to simulate from: ", Format (T,1,1), "\n");

characters = {{"A","C","G","T"}
              {"1", "", "", ""}};

DataSet testData = Simulate (T, eqf, characters, 1, 0);

for(site=0; site<100; site+=1){
  siteScale := siteRates[site__];
  DataSet sim = Simulate (T, eqf, characters, 1, 0);

  if (site == 0) {
    DataSet testData = Simulate (T, eqf, characters, 1, 0);
  } else {
    DataSet testData = Concatenate(testData, sim);
  }
}
DataSetFilter testDataF = CreateFilter(testData, 1);
fprintf(stdout, "\n", testDataF);
fprintf(stdout, "\n", siteRates);