LF_SMOOTHING_SCALER         = 0.1;

// namespace 'utility' for convenience functions 
LoadFunctionLibrary("lib2014/UtilityFunctions.bf");

// namespace 'io' for interactive/datamonkey i/o functions
LoadFunctionLibrary("lib2014/IOFunctions.bf");
LoadFunctionLibrary("lib2014/models/codon/MG_REV.bf");

// namespace 'estimators' for various estimator related functions
LoadFunctionLibrary("lib2014/tasks/estimators.bf");

io.displayAnalysisBanner ({"info" : "Estimate dS and dN from an alignment 
                            (pooling over branches and sites)",
                           "version" : "1.00",
                           "reference" : "https://github.com/veg/hyphy/issues/229",
                           "authors" : "Sergei L Kosakovsky Pond",
                           "contact" : "spond@ucsd.edu",
                           "requirements" : "in-frame codon alignment and a phylogenetic tree"         
                          } );

codon_data_info = utility.promptForGeneticCodeAndAlignment ("codon_data", "codon_filter");

LoadFunctionLibrary ("CodonTools.def");

tree 	  = utility.loadAnnotatedTopology (1);

DataSetFilter dnds_estimator.codon_data = CreateFilter (codon_filter, 3, "", "", codon_data_info["stop"]);

dnds_estimator.model =  model.generic.define_model ("models.codon.MG_REV.modelDescription", 
                                                         "dnds_estimator.mg_model", 
                                                         {"0" : parameters.quote("local"), "1" : codon_data_info["code"]}, 
                                                         "dnds_estimator.codon_data", 
                                                         None);

non_stop = (codon_data_info["code"])["_MATRIX_ELEMENT_VALUE_!=10"];
syn_site_counts = 3/(((_S_NS_POSITIONS_)[0][-1])[non_stop]*dnds_estimator.model["EFV"])[0];
nonsyn_site_counts = 3/(((_S_NS_POSITIONS_)[1][-1])[non_stop]*dnds_estimator.model["EFV"])[0];
                                                                                              
model.applyModelToTree ("dnds.tree", tree, {"default" : dnds_estimator.model} , None);

LikelihoodFunction dnds_estimator.likelihoodFunction = (dnds_estimator.codon_data, dnds.tree);
Optimize (dnds_estimator.mles, dnds_estimator.likelihoodFunction);

LoadFunctionLibrary ("dSdNTreeTools");

dnds_estimator.lengths = ReturnVectorsOfCodonLengths (ComputeScalingStencils(0), "dnds.tree");

dnds_estimator.lengths.syn    = +dnds_estimator.lengths["Syn"] * syn_site_counts;
dnds_estimator.lengths.nonsyn = +dnds_estimator.lengths["NonSyn"] * nonsyn_site_counts;

fprintf (stdout, "\ndS = ", dnds_estimator.lengths.syn, 
                 "\ndN = ", dnds_estimator.lengths.nonsyn,
                 "\n");
                 
tt = ReturnVectorsOfCodonLengths (ComputeScalingStencils(0), "dnds.tree");

keys = {{"Syn","NonSyn"}};

flatTree = dnds.tree ^ 0;

LoadFunctionLibrary ("TreeTools");

for (k = 0; k < Columns (keys); k+=1) {
    distanceAnnotator = {};
    for (b = 0; b < Abs (flatTree) - 1; b+=1) {
        distanceAnnotator [ (flatTree[b])["Name"] ] = (tt[keys[k]])[b] * dnds_estimator.codon_data.sites * 3;
    }   
    fprintf (stdout, "\n", keys[k], "\n", PostOrderAVL2StringDistances (flatTree, distanceAnnotator), "\n");
}