/* 
 * pipeline input parameters 
 */
params.reads = "$projectDir/data/test_data/*fastq.gz"
// params.reads = "$projectDir/data/test_data/data39/CDAN_1.fastq.gz"
params.db = "$projectDir/data/hla/fasta/{A,B,C}_gen.fasta"
params.db_class2 = "$projectDir/data/hla/fasta/D*_nuc.fasta"
params.primers = "$projectDir/data/primers/primers2.csv"
params.outdir = "$projectDir/results/"
params.g_groups = "$projectDir/data/hla_nom_g.txt"
params.threshold = 0
params.out = 'aggregated_results'

log.info """\
         HLATYPING - N F   P I P E L I N E    
         ===================================
         reads        : ${params.reads}
         primers      : ${params.primers}
         outfile      : ${params.out}
         outdir       : ${params.outdir}
         """
         .stripIndent()


// Process to index genome data
process Index {
    tag { genome.baseName }
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
    demultiplex.py $reads $params.primers "${reads.baseName.replace('.fastq', '')}_demultiplexed_output"
    """
}

// Function to extract gene type from a filename
String getGeneTypeFromName(String name) {
    String geneType = name.split(/[-_]/)[0]
    return geneType
}
// Function to extract gene type from a filename
String getGeneAndSampleFromName(String name) {
    String geneAndSample = name.split(/[-_]/)[0] + '_' + name.split(/[-_]/)[1]
    return geneAndSample
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
    kma -i $demultiplexed/${primer} -tmp -t_db $indexed -eq 10 -5p 20 -3p 20 -1t1 -bcNano -lc -proxi -0.98 -o '${primer.replace('.fastq', '_')}${demultiplexed.baseName.replace('_demultiplexed_output', '_')}types'
    """
}

// Process for sorting typing results
process SortTypingResults {
    tag { locus }
    input:
    tuple val(locus), path(typing_results)

    output:
    path "${typing_results.baseName}_sorted.txt", emit: sorted_results

    script:
    """
    last_col=\$(head -n 1 ${typing_results} | awk '{print NF}')
    cut -f 1,2 ${typing_results} | sort -k 5 -n -r | head -n 3 > ${typing_results.baseName}_sorted.txt
    """
}

// Process for comparing proportions between the second most deep read and positive and negative controls
process CompareProportions {
    errorStrategy = "ignore"
    tag { file }

    input:
    path file

    output:
    path "${file.baseName}_selected.txt", emit: selected_records

    script:
    """

    echo "Processing file: ${file}"
    awk 'NR <= 3 { score[NR] = \$NF; line[NR] = \$0 }
    END {
        diff12 = score[1]
        diff13 = score[2]
        if (diff13 / diff12 < ${params.threshold}) 
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
    tag { file }

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
    tag { geneType }

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
    kma -i ${geneType}_${fileName}_matched_sequences.fasta -tmp -t_db  $indexed -proxi -0.98 -o '${geneType}_${fileName}_results'
    """
}

// Process for collecting all results (for each locus) for one file (usually person)
process AggregateResults {
    tag { fileName }

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

process TranslateAllelesToGGroups {
    tag { fileName }

    input:
    tuple val(fileName), path(input_file)

    output:
    path "${fileName}_translated.txt"

    script:
    """
    translate_to_g.py --input_file $input_file --g_group_file $params.g_groups --output_file ${fileName}_translated.txt
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
    path params.out + '.csv'

    script:
    // call a transform_and_aggregate script with all csv files and a resulting file name
    """
    transform_and_aggregate.py ${csv_files} ${params.out}.csv
    """
}



workflow{
    // INDEXING
    indexed_ch = Channel.fromPath(params.db)
        .concat(Channel.fromPath(params.db_class2))
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
    .map{ it -> tuple(getGeneAndSampleFromName(it[1].baseName), it[1], it[2]) }
    // mapped = typing_res.map { typed -> tuple(typed[0], typed[1]) }
    // | view()

    // EXTRACTING 3 BEST MATCHES FOR EACH LOCUS, DECIDING IF HOMO OR HETEROZYGOUS AND EXTRACTING THE MATCHED SEQUENCES
    sorted_typing_res = SortTypingResults(typing_res.map{it[0..1]})
    typing_res_from_comparison = CompareProportions(sorted_typing_res)
    // | view()
    selected_records = ExtractNames(typing_res_from_comparison)
    .map{ file -> tuple(getGeneAndSampleFromName(file.baseName), file) }
    // | view()

    extraction_input = selected_records
    .join(typing_res)
    .map{it.getAt([0,1,3])}
    // | view()

    // indexed_hla_path = indexed_ch
    // .filter { it[0] == 'hla' }
    // .map { it[1] } // Maps the channel to get only the indexPath

    // EXTRACTING SEQUENCES MATCHING THE 1 OR 2 HIGHEST HITS
    matching_sequences_ch = ExtractMatchingSequences(extraction_input) 
    .map{ file -> tuple(getGeneTypeFromName(file.baseName), file) }
    .combine(indexed_ch, by:0) // Combining with the HLA index path
    // | view()

    // ensuring to have both - names of alleles and files
    file_names_ch = matching_sequences_ch.map{ it -> 
        def splitted_name = it[1].baseName.toString().split('_', 2)
        return [splitted_name[0], splitted_name[1].replace('_matched_sequences.fasta', '')]
        } 

    // FINAL MAPPING INPUT - SEQUENCES, ALLELE NAME, FILE NAME
    matching_sequences = matching_sequences_ch.join(file_names_ch)
    // | view()

    // MAPPING THE BEST RESULTS FROM KMA WITH THE WHOLE DATABASE
    final_mapping_res = FinalMapping(matching_sequences)
    .map { it ->
        def splitted_name = it.baseName.toString().split('_', 2)
        return [splitted_name[1].replace('_results.res', ''), it]
    }
    .groupTuple()
    // | view()

    // COLLECTING RESULTS FOR ALL LOCI
    aggregated_res_ch = AggregateResults(final_mapping_res)
    .map { it ->
        return [it.baseName.toString().replace('_matched_sequences_results_aggregated', ''), it]
    }
    // | view()

    translated_res_ch = TranslateAllelesToGGroups(aggregated_res_ch)
    // | view()

    // CONVERTING THE RESULTS TO A CSV FILE
    csv_files_ch = ConvertToCsv(translated_res_ch)
    // | view()

    // COLLECTING ALL CSV FILES TO CREATE A SINGLE CSV WITH ALL THE RESULTS
    TransformAndAggregateCSVs(csv_files_ch.collect())
    // | view()
}