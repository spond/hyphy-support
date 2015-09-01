ACCEPT_ROOTED_TREES = 1;

/*
DataSet       binary_data      = ReadDataFile("data/binary.seq");
DataSetFilter binary_filter    = CreateFilter(raw_data,1);
*/


general_binary_rate_matrix     = 
    {{*, r01}
     {r10, *}};
     
global freq_0 = 0.8;
freq_0        :< 1;

root_frequencies = {{freq_0}{1-freq_0}};

Model general_binary_model = (general_binary_rate_matrix, root_frequencies, 0);

binary_characters = {{"0","1"}{"1",""}};

Tree binary_tree = (317,((((135,(529,105r)),(719,136)),6760,((113,9939),(256,(822,159)))),6767));

bnames = BranchName (binary_tree, -1);

for (k = 0; k < Columns (bnames) - 1; k+=1) {
    this_branch = bnames[k];
    Eval ("binary_tree.`this_branch`.r01 = 0.2");
    Eval ("binary_tree.`this_branch`.r10 = 0.8");
}

DataSet 		sim  = Simulate     (binary_tree,root_frequencies,binary_characters,1000,0);
DataSetFilter   sim_f = CreateFilter (sim, 1);

DATA_FILE_PRINT_FORMAT = 9;
DATAFILE_TREE = Format (binary_tree,0,0);
IS_TREE_PRESENT_IN_DATA = 1;

fprintf ("data/binary.seq", CLEAR_FILE, sim_f);

/*
LikelihoodFunction LikFn = (filt_data, Tree01);
Optimize (paramValues, LikFn);
fprintf (stdout, LikFn);
*/