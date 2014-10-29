ACCEPT_ROOTED_TREES         = 1;
LF_SMOOTHING_SCALER         = 0.1;


DataSet       binary_data      = ReadDataFile("data/binary.seq");
DataSetFilter binary_filter    = CreateFilter(binary_data,1);

/* define a rate matrix with 'local', i.e. branch specific rates
   parameters (0-1 and 1-0) */
   
general_binary_rate_matrix     = 
    {{*, r01}
     {r10, *}};
   
/* define the 'global', i.e. alignment wide frequencies at the root */
     
global freq_0 = 0.3;
freq_0        :< 1;
root_frequencies = {{freq_0}{1-freq_0}};

/* define the model */

Model general_binary_model = (general_binary_rate_matrix, root_frequencies, 0);

/* set up and optimize the likelihood function */

Tree binary_tree = DATAFILE_TREE;
LikelihoodFunction LikFn = (binary_filter, binary_tree);
Optimize (mles, LikFn);


fprintf (stdout, "
    Log (L) = ", mles[1][0], " 
    Inferred root frequencies
        f (0) : ", Format (freq_0, 5, 2),"
        f (1) : ", Format (1-freq_0, 5, 2), "\n");
        
/* this code is a hack to display estimated frequncies of the '0' character 
   in the phylogenetic tree */
   
LoadFunctionLibrary ("TreeTools");
flat_tree = binary_tree ^ 0;
flat_tree ["annotate"][""];

fprintf (stdout, "\nAnnotated tree ([freq 0 / freq 1])\n",
            PostOrderAVL2StringAnnotate (flat_tree, 1, "weights"), "\n");

function annotate (key, value) {
    if (value["Parent"]) {
        n_name = "binary_tree." + value["Name"];
        n_r01     = Eval ("`n_name`.r01");
        n_r10     = Eval ("`n_name`.r10");
        
        if (n_r01 + n_r10) {
            f0 = n_r01 / (n_r01 + n_r10);
            value ["weights"] = "" + Format (f0,0,3) + "/" + Format(1-f0,0,3);
            value ["Length"]  = f0 * n_r01 + (1-f0) * n_r01;
        } else {
            value ["weights"] = "Undefined";
        }
        
        
    }
}