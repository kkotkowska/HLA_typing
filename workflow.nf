/* 
 * pipeline input parameters 
 */
params.reads = "$projectDir/data/test_data/*fastq.gz"
// params.reads = "$projectDir/data/test_data/data39/CDAN_1.fastq.gz"
params.db = "$projectDir/data/hla/fasta/*_gen.fasta"
params.primers = "$projectDir/data/primers/primers1.csv"
params.outdir = "$projectDir/results/"

log.info """\
         HLATYPING - N F   P I P E L I N E    
         ===================================
         genome       : ${params.db}
         reads        : ${params.reads}
         primers      : ${params.primers}
         out          : ${params.outdir}
         """
         .stripIndent()

 
// Process to demultiplex FASTQ files based on primer sequences
process Demultiplex {
    tag {reads.baseName}

    input:
    path reads

    output:
    path "${reads.baseName.replace('.fastq', '')}_demultiplexed_output", emit: demultiplexed_dir

    script:
    // Call demultiplex.py script with the reads, primers, and output directory
    """
    echo $params.primers
    demultiplex.py $reads $params.primers "${reads.baseName.replace('.fastq', '')}_demultiplexed_output"
    """
}

// Function to extract gene type from a filename
String getGeneTypeFromName(String name) {
    String geneType = name.split(/[-_]/)[0]
    return geneType
}

// Process to index genome data
process Index {
    input:
    path genome

    output:
    path "indexed/${genome.baseName}_indexed.*"

    script:
    // Create an indexed directory and perform indexing using kma
    """
    mkdir indexed
    kma index -i $genome -m 14 -o indexed/${genome.baseName}_indexed
    """
}

// Process for typing based on demultiplexed data and indexed database
process Typing {
    tag { primer }

    input:
    tuple val(primer), path(demultiplexed), val(indexed)
    
    output:
    tuple val({primer.replace('.fastq', '')}), path("${primer.replace('.fastq', '_')}${demultiplexed.baseName.replace('_demultiplexed_output', '_')}types.res"), path("${primer.replace('.fastq', '_')}${demultiplexed.baseName.replace('_demultiplexed_output', '_')}types.fsa")

    
    script:
    // Call kma for typing
    """
    kma -i $demultiplexed/${primer} -tmp -t_db $indexed -eq 10 -5p 20 -3p 20 -1t1 -bc -bcNano -lc -proxi -0.98 -o '${primer.replace('.fastq', '_')}${demultiplexed.baseName.replace('_demultiplexed_output', '_')}types'
    """
}

// Process for sorting typing results
process SortTypingResults {
    input:
    tuple val(locus), path(typing_results)

    output:
    path "${typing_results.baseName}_sorted.txt", emit: sorted_results

    script:
    """
    last_col=\$(head -n 1 ${typing_results} | awk '{print NF}')
    cut -f 1,9 ${typing_results} | sort -k 5 -n -r | head -n 3 > ${typing_results.baseName}_sorted.txt
    """
}

// Process for comparing proportions between the second most deep read and positive and negative controls
process CompareProportions {
    tag "Comparing proportions in ${file}"

    input:
    path file

    output:
    path "${file.baseName}_selected.txt", emit: selected_records

    script:
    """
    echo "Processing file: ${file}"
    awk 'NR <= 3 { depth[NR] = \$NF; line[NR] = \$0 }
    END {
        diff12 = depth[1] - depth[2]
        diff13 = depth[1] - depth[3]
        if (diff13 / diff12 < 1.1) 
            print line[1] > "${file.baseName}_selected.txt"
        else {
            print line[1] > "${file.baseName}_selected.txt"
            print line[2] >> "${file.baseName}_selected.txt"
        }
    }' ${file}
    """
}

// Process for extracting names for further analysis
process ExtractNames {
    tag "Extracting names from ${file}"

    input:
    path file

    output:
    path "${file.baseName}_names.txt", emit: names_list

    script:
    """
    awk -F '\\t' '{print \$1}' ${file} > "${file.baseName}_names.txt"
    """
}

// Process for extracting sequences from a resulting file from KMA only for the hits
process ExtractMatchingSequences {
    tag "${geneType}"

    input:
    tuple val(geneType), path(namesFile), path(fsaFile)

    output:
    path("${namesFile.baseName.replace('_types_sorted_selected_names', '')}_matched_sequences.fasta")

    script:
    // call extract_typed_sequences script with names, sequences, and a resulting file name
    """
    extract_typed_sequences.py $namesFile $fsaFile ${namesFile.baseName.replace('_types_sorted_selected_names', '')}_matched_sequences.fasta
    """
}

// Process for typing of the extracted sequences
process FinalMapping {
    tag "${geneType}"

    input:
    tuple val(geneType), path("${geneType}_${fileName}_matched_sequences.fasta"), val(indexed), val(fileName)

    output:
    path("${geneType}_${fileName}_results.res")

    script:
    """
    kma -i ${geneType}_${fileName}_matched_sequences.fasta -tmp -t_db  $indexed -5p 20 -3p 20 -bc -bcNano -proxi -0.98 -1t1 -o '${geneType}_${fileName}_results'
    """
}

// Process for collecting all results (for each locus) for one file (usually person)
process AggregateResults {
    input:
    tuple val(fileName), path(files)

    output:
    path("${fileName}_aggregated.txt")

    script:
    // The script concatenates the contents of all files, keeping the header from the first file only.
    """
    HEADER=true
    for file in \$(echo ${files}); do
        cut -f 1 \$file >> ${fileName}_aggregated.txt
    done
    """
}

// Process for converting the results to a csv
process ConvertToCsv {
    input:
    path input_file

    output:
    path "${input_file.baseName}.csv"

    script:
    // call get_results_csv script with an input file and a resulting file name
    """
    get_results_csv.py ${input_file} ${input_file.baseName}.csv
    """
}

// Process for creating a csv file with all results (for all examined files)
process TransformAndAggregateCSVs {
    publishDir params.outdir

    input:
    path csv_files

    output:
    path "aggregated_results.csv"

    script:
    // call a transform_and_aggregate script with all csv files and a resulting file name
    """
    transform_and_aggregate.py ${csv_files} aggregated_results.csv
    """
}



workflow{
    // INDEXING
    indexed_ch = Channel.fromPath(params.db) 
        | Index 
        | map {it.get(2)}
        | map { indexed -> tuple(getGeneTypeFromName(indexed.baseName), indexed.toString().replace('.name','')) }
        // | view()

    reads_ch = Channel
    .fromPath(params.reads)
    .ifEmpty { error "No reads found!" }

    // DEMULTIPLEXING
    demultiplexed = Demultiplex(reads_ch)
    demultiplexed_files_ch = demultiplexed.flatMap { dir ->
        dir.list().collect { file ->
        tuple(file.replace('.fastq', ''), file, dir)
        }
    }

    // TYPING INPUT - DEMULTIPLEXED + INDEXED DB
    typing_input = demultiplexed_files_ch.combine(indexed_ch, by:0)
    .map{ key, fastq, dir, db -> [fastq, dir, db] }
    // | view()

    // TYPING
    typing_res = Typing(typing_input)
    mapped = typing_res.map { typed -> tuple(typed[0], typed[1]) }
    // | view()

    // EXTRACTING 3 BEST MATCHES FOR EACH LOCUS AND COMBINING WITH THE FULL INDEXED DATABASE 
    sorted_typing_res = SortTypingResults(mapped)
    typing_res_from_comparison = CompareProportions(sorted_typing_res)
    selected_records = ExtractNames(typing_res_from_comparison)
    .map{ file -> tuple(getGeneTypeFromName(file.baseName), file) }
    // | view()

    extraction_input = selected_records
    .join(typing_res)
    .map{it.getAt([0,1,3])}

    indexed_hla_path = indexed_ch
    .filter { it[0] == 'hla' }
    .map { it[1] } // Maps the channel to get only the indexPath

    // EXTRACTING SEQUENCES MATCHING THE 1 OR 2 HIGHEST HITS
    matching_sequences_ch = ExtractMatchingSequences(extraction_input) 
    .map{ file -> tuple(getGeneTypeFromName(file.baseName), file) }
    .combine(indexed_hla_path) // Combining with the HLA index path
    // | view()

    // ensuring to have both - names of alleles and files
    file_names_ch = matching_sequences_ch.map{ it -> 
        def splitted_name = it[1].baseName.toString().split('_', 2)
        return [splitted_name[0], splitted_name[1].replace('_matched_sequences.fasta', '')]
        } 
    // | view()

    // FINAL MAPPING INPUT - SEQUENCES, ALLELE NAME, FILE NAME
    matching_sequences = matching_sequences_ch.join(file_names_ch)
    // | view()

    // MAPPING THE BEST RESULTS FROM KMA WITH THE WHOLE DATABASE
    final_mapping_res = FinalMapping(matching_sequences)
    // | view()
    .map { it ->
        def splitted_name = it.baseName.toString().split('_', 2)
        return [splitted_name[1].replace('_results.res', ''), it]
    }
    .groupTuple()
    // | view()

    // COLLECTING RESULTS FOR ALL LOCI
    aggregated_res_ch = AggregateResults(final_mapping_res)
    // | view()

    // CONVERTING THE RESULTS TO A CSV FILE
    csv_files_ch = ConvertToCsv(aggregated_res_ch)

    // COLLECTING ALL CSV FILES TO CREATE A SINGLE CSV WITH ALL THE RESULTS
    TransformAndAggregateCSVs(csv_files_ch.collect())
    | view()
}