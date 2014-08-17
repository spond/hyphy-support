vector = {{10,4,-5,21,0}};

fprintf (stdout, "\nVector min = ", 
                 Min (vector, 0), 
                 "\nVector max = ", 
                 Max (vector, 0), 
                 "\n");
             
fprintf (stdout, "\nVector max with index = ", Max (vector, 1), "\n");


LoadFunctionLibrary ("ReadDelimitedFiles");

csv_data = ReadCSVTable ("test.csv", 1);

fprintf (stdout, "\nColumn headers:", csv_data[0],
                 "\nTable data:", csv_data[1], "\n");
                 
                 
LoadFunctionLibrary ("ProbabilityDistributions");

// populate a vector with standard normals
// the syntax below simply applies the formula 'sampleFromNormal'
// to each entry in the matrix

random_normal_values = {1000, 1} ["sampleFromNormal()"]; 

mu = 10;
sigma = 0.1;

// transform the deviates using matrix arithmetic

random_normal_values = random_normal_values * sigma + mu;

// confirm that we got what we wanted

LoadFunctionLibrary ("DescriptiveStatistics");

PrintDescriptiveStats ("Random deviates with mean " + mu + " and std.dev " + sigma,
                        GatherDescriptiveStats (random_normal_values));
                        
         
// lfunction will give the variables inside the function a local scope         
                        
lfunction GaussianDensity (x, mu, sigma) {
    return Exp (- (x-mu)^2 * 0.5 / sigma^2) / sigma / Sqrt (8*Arctan (1));
}

function GaussianPercentile (prob, mu, sigma) {
    // ZCDF is the standard normal cumulative distribution function
    // FindRoot is documented here http://hyphy.org/w/index.php/FindRoot
    FindRoot (z, ZCDF ((t-mu)/sigma) - prob, t, -100, 100);
    return z;
     
}

for (x = 0; x < 5; x+=0.5) {
    fprintf (stdout, "Normal (", x, ", mu = ", 2, ", sigma = ", 1.5, ") = ",
                      GaussianDensity (x, 2, 1.5), "\n");
}

for (x = 0.05; x < 1; x+=0.05) {
    fprintf (stdout, Format (x*100, 5, 0), " percentile of Normal (mu = ", 2, ", sigma = ", 1.5, ") = ",
                      GaussianPercentile (x, 2, 1.5), "\n");
}