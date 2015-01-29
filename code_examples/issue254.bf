function define_logistic_mapper (dimension, model_name, logistic_regression_parameter_prefix) {
    // see http://en.wikipedia.org/wiki/Logistic_regression
    
    // zero NxN matrix   
    define_logistic_mapper.matrix = model_name + ".transition_matrix";    
    ExecuteCommands ("`define_logistic_mapper.matrix` = {dimension, dimension};");

    // fake Nx1 frequencies matrix [needed to define the model, but will be ignored]   
    define_logistic_mapper.frequencies = model_name + ".frequencies";
    ExecuteCommands ("`define_logistic_mapper.frequencies` = {dimension,1}[\"1/dimension\"]"); // set all entries to 1/dimension
    
    // logistic model parameters    
    define_logistic_mapper.intercept = logistic_regression_parameter_prefix + ".beta0";
    define_logistic_mapper.slope     = logistic_regression_parameter_prefix + ".beta1";
    
    ExecuteCommands ("global `define_logistic_mapper.intercept` = 0.0");
    *define_logistic_mapper.intercept :> -1e10;
    *define_logistic_mapper.intercept :< 1e10;

    ExecuteCommands ("global `define_logistic_mapper.slope` = 0.1");
    *define_logistic_mapper.slope :> -1e10;
    *define_logistic_mapper.slope :< 1e10;
    
    for ( define_logistic_mapper.index = 0;  define_logistic_mapper.index < dimension; define_logistic_mapper.index += 1) {
    
        ExecuteCommands ("`define_logistic_mapper.matrix`[define_logistic_mapper.index][define_logistic_mapper.index] = 0");
        
        ExecuteCommands ("`define_logistic_mapper.matrix`[define_logistic_mapper.index][1] := " + 
                "1/(1+Exp (-`define_logistic_mapper.intercept` - `define_logistic_mapper.slope` * (" + (define_logistic_mapper.index - dimension$2) + ")))") ;
        
        // probability of positive phenotype given "hidden value = index"
 
        ExecuteCommands ("`define_logistic_mapper.matrix`[define_logistic_mapper.index][0] := " + 
                "1- 1/(1+Exp (-`define_logistic_mapper.intercept` - `define_logistic_mapper.slope` * (" + (define_logistic_mapper.index - dimension$2) + ")))") ;
            
        // 1 - probability of positive phenotype given "hidden value = index"
    }
    
    ExecuteCommands ("Model `model_name` = (\"`define_logistic_mapper.matrix`\", `define_logistic_mapper.frequencies`, EXPLICIT_FORM_MATRIX_EXPONENTIAL)");
    
}

function tree_extender (tree_id, leaf_suffix, model_name) {
    tree_extender.leaf_count = TipCount (*tree_id);
    for (tree_extender.counter = 0; tree_extender.counter < tree_extender.leaf_count; tree_extender.counter += 1) {
        tree_extender.leaf_name = TipName (*tree_id, tree_extender.counter);
        
        
        (*tree_id) + {"WHERE": tree_extender.leaf_name,
                      "NAME": tree_extender.leaf_name + leaf_suffix};
                      
        ExecuteCommands ("SetParameter (" + tree_id + "." + tree_extender.leaf_name + leaf_suffix + ", MODEL, `model_name`)");
    }
}


// AN EXAMPLE OF USE

define_logistic_mapper (128, "LOGISTIC_MAPPER", "LOGISTIC_MAPPER");

INCLUDE_MODEL_SPECS = 1;
UseModel (USE_NO_MODEL);
Tree T = ((a,b),(c,d),e);

fprintf (stdout, "BEFORE\n", T, "\n");

tree_extender ("T", "_phenotype", "LOGISTIC_MAPPER");

fprintf (stdout, "AFTER\n", T, "\n");