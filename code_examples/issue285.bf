// namespace 'utility' for convenience functions 
LoadFunctionLibrary("lib2014/UtilityFunctions.bf");

// namespace 'io' for interactive/datamonkey i/o functions
LoadFunctionLibrary("lib2014/IOFunctions.bf");
LoadFunctionLibrary("lib2014/models/DNA/GTR.bf");

// namespace 'estimators' for various estimator related functions
LoadFunctionLibrary("lib2014/tasks/estimators.bf");

io.displayAnalysisBanner ({"info" : "Fit a shared GTR + Beta/Gamma model to multiple partitions 
                            (different sizes/trees)",
                           "version" : "1.00",
                           "reference" : "https://github.com/veg/hyphy/issues/285",
                           "authors" : "Sergei L Kosakovsky Pond",
                           "contact" : "spond@ucsd.edu",
                           "requirements" : "one or more nucleotide alignments with trees"         
                          } );

list_of_files = io.get_a_list_of_files ("data/file.lst"); 
// prompt for a list from files from the console or
// do io.get_a_list_of_files (file_name_with_list ); 

io.reportProgressMessage ("Multi-file", "Read " + Abs (list_of_files) + " file paths");

// now read the input files and trees

input_files  = {};
tree_strings = {};
filter_set   = {Abs (list_of_files), 1};

function read_file_callback (id, file_path) {
    input_files [id] = utility.readNucleotideAlignment (file_path, "data_" + id, "filter_" + id);
    io.reportProgressMessage ("Multi-file", "Loaded " + (input_files [id])["sequences"] + " sequences with " + (input_files [id])["sites"] + " sites from `file_path`");
    io.reportProgressMessage ("Multi-file", "Loading tree for file `file_path`");
    filter_set  [0+id] = (input_files [id])["datafilter"];
    tree_strings[id] = io.getTreeString (1);
}


list_of_files ["read_file_callback"][""];

fprintf (stdout, filter_set, "\n");

// define a combined GTR model

shared.model =  model.generic.define_model ("models.DNA.GTR.modelDescription", // this is the model template [function ] defined in lib2014/models/DNA/GTR.bf and passed by reference
                                            "shared.GTR", // this is the namespace for Q matrices, global parameters, etc
                                            {"0" : "terms.global"}, // this is the set of arguments to be passed to the template definition
                                            filter_set, // this is the set of data filters to base empirical frequencies on 
                                            None // if not None, will override the default (empirical in this case)
                                            ); 

fprintf (stdout, shared.model, "\n");

return;

model.applyModelToTree ("estimators.fitGTR.tree", tree, {"default" : estimators.fitGTR.model} , None);

LikelihoodFunction estimators.fitGTR.likelihoodFunction = (estimators.fitGTR.nuc_data, estimators.fitGTR.tree);

if (Type (initial_values) == "AssociativeList") {
    utility.toggleEnvVariable ("USE_LAST_RESULTS", 1);
    estimators.applyExistingEstimates ("estimators.fitGTR.likelihoodFunction", {"estimators.fitGTR.gtr" : estimators.fitGTR.model}, initial_values);
}

Optimize (estimators.fitGTR.mles, estimators.fitGTR.likelihoodFunction);
if (Type (initial_values) == "AssociativeList") {
    utility.toggleEnvVariable ("USE_LAST_RESULTS", None);
}
    
estimators.fitGTR.results = estimators.extractMLEs ("estimators.fitGTR.likelihoodFunction", {"estimators.fitGTR.gtr" : estimators.fitGTR.model});

estimators.fitGTR.results["LogL"]            = estimators.fitGTR.mles[1][0];
estimators.fitGTR.results["parameters"] = estimators.fitGTR.mles[1][1] + 3;

DeleteObject (estimators.fitGTR.likelihoodFunction);

return estimators.fitGTR.results;

