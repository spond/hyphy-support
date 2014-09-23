SetDialogPrompt ("Please specify a nucleotide or amino-acid data file:");

DataSet 	  msa 		   = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData   = CreateFilter (msa,1);

SelectTemplateModel(filteredData);

ACCEPT_ROOTED_TREES = 1;
Tree   user_tree = DATAFILE_TREE;

ChoiceList (reconstruction_type,
            "Root sequence reconstruction",1,SKIP_NONE,
            "Joint", "Ancestral characters are chosen to maximize the conditional probability over the entire tree",
            "Marginal", "Ancestral characters are chosen to maximize the probability at a given node (marginalizing over all other nodes)");
           
if (reconstruction_type < 0) {
    return;
}          
  
LikelihoodFunction lf = (filteredData, user_tree);
Optimize (lf.mles, lf);

if (reconstruction_type == 0) {
    DataSet	 				ancestors = ReconstructAncestors (lf);
} else {
    DataSet	 		        ancestors = ReconstructAncestors (lf,MARGINAL);
}

DataSetFilter ancestral_filter = CreateFilter (ancestors, 1);

// the root sequence is always the first one in the file

GetDataInfo (root_ancestor, ancestral_filter, 0);

fprintf (stdout, "Inferred the `SELECTION_STRINGS` root sequence as \n", root_ancestor, "\n");

