codons_to_replace = {"GGG" : 1};

DataSet         nuc_data_to_filter = ReadDataFile ("data/hiv.seq");
DataSetFilter   nuc_data_to_filterF = CreateFilter (nuc_data_to_filter, 1);

filtered_FASTA = ""; filtered_FASTA * (nuc_data_to_filter.species * nuc_data_to_filter.sites); 
// filtered_FASTA is now a string buffer, initial string allocation is the number of characters in the file we just read in

for (seq = 0; seq < nuc_data_to_filter.species; seq += 1) {
    GetString   (seq_name, nuc_data_to_filterF, seq);
    GetDataInfo (seq_data, nuc_data_to_filterF, seq);
    filtered_FASTA * (">`seq_name`\n"); // push the FASTA header for this sequence onto string buffer
    for (char = 0; char < nuc_data_to_filterF.sites; char += 3) {
        current_codon = seq_data[char][char+2];
        if (codons_to_replace[current_codon]) {
            filtered_FASTA * "NNN";
        } else {
            filtered_FASTA * current_codon;
        }   
    }
    filtered_FASTA * "\n";
}

filtered_FASTA * 0; // close the buffer

fprintf (stdout, "\n", filtered_FASTA, "\n");

DataSet codons_replaced = ReadFromString (filtered_FASTA);
DataSetFilter codons_replacedF = CreateFilter (codons_replaced, 1);

