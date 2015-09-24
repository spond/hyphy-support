function file_exists (filename) {
    return !filename;
}

// D = number of points; S = number of sites;
// G = number of grid samples set globally

function compute_sites_posteriors (siteLL, grid, grid_weights) {
   // pr (a,b | site) = pr (site | a,b) * pr (a,b) / pr (site)


    // for each site i, grid sample j, and grid point k we need to compute
    // pr (site_i | a = a^j_k, b = b^j_k) pr (a = a^j_k, b = b^j_k)
    // this is going to be a list (indexed on i) of D x G matrices
    // with element [k][m] being pr (site_k | a_m,b_m)

    _full_site_posteriors = {};

    for (site = 0; site < S; site += 1) {
        _site_posteriors = {D,G};
        _site_conditionals = siteLL [-1][site];
            // this extracts a Dx1 column of conditional site probs at all grid points

        _site_conditionals = Transpose(_site_conditionals * ({1,G} ["1"]));
            // make a DxG matrix where each column "i" holds G copies of site likelihood
            // conditioned on grid point i

        _unscaled_site_posteriors = _site_conditionals $ grid_weights;

        // grid_weight [i][j] is the weight of the j-th grid point in the i-th sample
        // _unscaled_site_posteriors[i][j] = Pr (site | grid_point_i) * grid_weight[i][j]

        _normalizing_coefficients = _unscaled_site_posteriors * ({D,1}["1"]);



        posterior_ab = (_unscaled_site_posteriors * grid) / (_normalizing_coefficients * ({1,2}["1"]));
                    //  G x D * D x 2 (G x 2)

        _full_site_posteriors + posterior_ab;
    }

    return _full_site_posteriors;


    posterior_scratch = grid_weights * siteLL;
        // this computes \sum_{grid_sample_j} pr (site_i | a,b) pr (a, b), i.e.
        // element (i,j) of the matrix holds \sum_grid pr (site | a,b) * pr (a,b)
        // for grid sample i at site j; i.e. pr (site_i) given grid weight sample j

    posterior = {S, 2};
        // the posterior mean of (alpha, beta) for this sample at each site


    //fprintf (stdout, Rows (posterior_scratch), "x", Columns (posterior_scratch), "\n");
}

LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("ReadDelimitedFiles");
LoadFunctionLibrary ("WriteDelimitedFiles");
LoadFunctionLibrary ("DescriptiveStatistics");

SetDialogPrompt ("Please select a .grid_info file:");
fscanf (PROMPT_FOR_FILE, "NMatrix,Raw", grid_points, site_likelihoods);

D = Rows (grid_points); // number of grid points

site_likelihoods = Eval (site_likelihoods);

S = Columns (site_likelihoods ["conditionals"]);

FUBAR_path = splitFilePath  (LAST_FILE_PATH);

FUBAR_path_trunk = FUBAR_path["DIRECTORY"] +
                   DIRECTORY_SEPARATOR +
                   FUBAR_path["FILENAME"];

chain_count = 0;
grid_samples = {};
test_file_name = (FUBAR_path_trunk + ".samples." + chain_count);

while (file_exists(test_file_name)) {
    fscanf (test_file_name, "NMatrix,NMatrix", ignore, grid_sample);
    grid_samples + grid_sample;
    chain_count += 1;
    test_file_name = (FUBAR_path_trunk + ".samples." + chain_count);
}

assert (chain_count > 0, "Could not find any .samples files");

fprintf (stdout, "\n[LOADED SAMPLES FROM ", chain_count, " CHAINS]\n");

G = Rows (grid_samples[0]);

if (chain_count > 1) {
    raw_G = G;
    G = G $ chain_count * chain_count;
    fprintf (stdout, "\n[THINNING CHAIN SAMPLES DOWN TO ", G, "]\n");
    final_sample = {G, D};
    thinnned_index = 0;
    for (chain_id = 0; chain_id < chain_count; chain_id += 1) {
        for (sample_id = 0; sample_id < raw_G; sample_id += 1) {
            if (sample_id % chain_count == 0) {
                for (grid_point = 0; grid_point < D; grid_point += 1) {
                    final_sample [thinnned_index][grid_point] = (grid_samples[chain_id])[sample_id][grid_point];
                }
                thinnned_index += 1;
            }
        }
    }

} else {
    final_sample = grid_samples[0];
}

posteriors = compute_sites_posteriors ( site_likelihoods["conditionals"],
                           grid_points,
                           final_sample);

summary_statistics = {S, 5}; // index, mean alpha, var alpha, mean beta, var beta

for (s = 0; s < S; s += 1) {
    summary_statistics[s][0] = s+1;
    alpha_stats = GatherDescriptiveStats ((posteriors [s])[-1][0]);
        // there are many other stats available here: median, 2.5% etc
    summary_statistics[s][1] = alpha_stats["Mean"];
    summary_statistics[s][2] = alpha_stats["Variance"];

    beta_stats = GatherDescriptiveStats ((posteriors [s])[-1][1]);
        // there are many other stats available here: median, 2.5% etc
    summary_statistics[s][3] = beta_stats["Mean"];
    summary_statistics[s][4] = beta_stats["Variance"];
}

SetDialogPrompt ("Write summary .csv to:");
WriteSeparatedTable ("", {{"Codon","E[alpha]", "Var[alpha]", "E[beta]", "Var[beta]"}}, summary_statistics, 0, ",");
