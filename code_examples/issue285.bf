// namespace 'utility' for convenience functions 
LoadFunctionLibrary("lib2014/UtilityFunctions.bf");

// namespace 'io' for interactive/datamonkey i/o functions
LoadFunctionLibrary("lib2014/IOFunctions.bf");
LoadFunctionLibrary("lib2014/models/DNA/GTR.bf");

// namespace 'estimators' for various estimator related functions
LoadFunctionLibrary("lib2014/tasks/estimators.bf");

// namespace 'rate_variation' for various rate variation related functions
LoadFunctionLibrary("lib2014/models/rate_variation.bf");

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
filter_set   = {Abs (list_of_files), 1}; // filter id

function read_file_callback (id, file_path) {
    input_files [id] = utility.readNucleotideAlignment (file_path, "data_" + id, "filter_" + id);
    io.reportProgressMessage ("Multi-file", "Loaded " + (input_files [id])["sequences"] + " sequences with " + (input_files [id])["sites"] + " sites from `file_path`");
    io.reportProgressMessage ("Multi-file", "Loading tree for file `file_path`");
    filter_set  [0+id] = (input_files [id])["datafilter"];
    tree_strings[id] = utility.loadAnnotatedTopology (1);
}

list_of_files ["read_file_callback"][""];


// define a combined GTR model with rate variation

function gtr_i_gamma (options) {
    return rate_variation.add (models.DNA.GTR.modelDescription (terms.global),
                              options);
}

shared.model =  model.generic.define_model ("gtr_i_gamma", // this is the model template [function ] defined above
                "shared.GTR", // this is the namespace for Q matrices, global parameters, etc
                {"0" : 
                    {"type" : "Gamma+I",
                    "bins" : "4"}
                }, // this is the set of arguments to be passed to the template definition
                filter_set, // this is the set of data filters to base empirical frequencies on 
                None // if not None, will override the default (empirical in this case)
                );


// define a separate tree for each partition

tree_list = {Abs (list_of_files), 1};
function apply_model_callback (id, file_path) {
    shared.tree.id = "shared.tree." + id;
    model.applyModelToTree (shared.tree.id, tree_strings[0+id] ,  {"default" : shared.model}, None);
    tree_list [0+id] = shared.tree.id;
}   
list_of_files ["apply_model_callback"][""];

shared.results = estimators.fitLF (filter_set,
                                              tree_list,
                                              {"shared.GTR": shared.model} , // model ID -> model matrix map, needed to do branch evaluations etc
                                              None // no initial conditions
                                              );

utility.json_spool (shared.results, None);