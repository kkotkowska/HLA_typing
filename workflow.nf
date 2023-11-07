/* 
 * pipeline input parameters 
 */
params.reads = "$projectDir/data/test_data/Test1.fastq.gz"
// params.reads = "$projectDir/data/test_data/data39/CDAN_1.fastq.gz"
params.db = "$projectDir/data/imgt-hla/fasta/*_gen.fasta"
params.db_class2 = "$projectDir/data/imgt-hla/D*_nuc.fasta"
params.exon_file = "$projectDir/data/exons/*.fa"
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

 
process Demultiplex {
    input:
    path reads
    path primers

    output:
    path 'demultiplexed_output', emit: demultiplexed_dir

    script:
    """
    demultiplex.py $reads $primers
    """
}

String getGeneTypeFromName(String name) {
    // Split the name on the '-' and take the first element

    String geneType = name.split(/[-_]/)[0]
    return geneType
}

process Index {
    
    input:
    path genome

    output:
    path "indexed/${genome.baseName}_indexed.*"

    script:
    """
    mkdir indexed
    kma index -i $genome -m 14 -o indexed/${genome.baseName}_indexed
    """
}

process Typing {
    tag { primer }

    input:
    tuple val(primer), val(indexed), path(demultiplexed), path(outdir)
    
    output:
    tuple val({primer.replace('.fastq', '')}),path("${primer.replace('.fastq', '')}_types.res"), path("${primer.replace('.fastq', '')}_types.fsa")

    
    script:
    """
    kma -i $demultiplexed/${primer} -tmp -t_db $indexed -eq 10 -5p 20 -3p 20 -1t1 -bc -bcNano -lc -o '${primer.replace('.fastq', '')}_types'
    """
}

process SortTypingResults {

    input:
    tuple val(locus), path(typing_results)

    output:
    path "${locus}_sorted.txt", emit: sorted_results

    script:
    """
    last_col=\$(head -n 1 ${typing_results} | awk '{print NF}')
    cut -f 1,9 ${typing_results} | sort -k 5 -n -r | head -n 3 > ${locus}_sorted.txt
    """
}

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

process ExtractMatchingSequences {
    tag "${geneType}"

    input:
    tuple val(geneType), path(namesFile), path(fsaFile)

    output:
    path("${geneType}_matched_sequences.fasta")

    script:
    """
    extract_typed_sequences.py $namesFile $fsaFile ${geneType}_matched_sequences.fasta

    """
}

process FinalMapping {
    tag "${geneType}"

    input:
    tuple val(geneType), path("${geneType}_matched_sequences.fasta"), val(indexed)

    output:
    path("${geneType}_typing_results.res")

    script:
    """
    kma -i ${geneType}_matched_sequences.fasta -tmp -t_db $indexed  -eq 10 -5p 20 -3p 20 -1t1 -o '${geneType}_typing_results'
    """
}


process GatherResults {
    publishDir "results/"
    input:
    path results_list

    output:
    path "final_results.txt"

    script:
    """
    cut -f 1 ${results_list}  > final_results.txt
    """
}
process ConvertToCsv {
    tag "Converting to CSV"
    publishDir params.outdir

    input:
    path input_file

    output:
    path "${input_file.baseName}.csv"

    script:
    """
    get_results_csv.py ${input_file} ${input_file.baseName}.csv
    """
}


workflow{
    // DEMULTIPLEXING
    demultiplexed = Demultiplex(channel.of(params.reads), channel.of(params.primers))
    demultiplexed_files_ch = demultiplexed.flatMap { dir -> tuple(dir.list()) }  // List the contents of the directory
    .map { file -> tuple((file.replace('.fastq', '')), file)}
    // | view()

    // INDEXING
    indexed_ch = Channel.fromPath(params.db) 
        | Index 
        | map {it.get(2)}
        | map { indexed -> tuple(getGeneTypeFromName(indexed.baseName), indexed.toString().replace('.name','')) }
        // | view()

    // COMBINING BINS AND DATASET - TYPING INPUT
    concatenated = demultiplexed_files_ch.concat(indexed_ch)
               .groupTuple()
               .flatMap { key, values -> 
                  def tuples = []

                  // Only process keys with matched values
                  if(values.size() > 1) {
                      values.each { value ->
                          tuples.add([key, value])
                      }
                  }
                  return tuples
               }
    // Splitting data
    fastq_ch = concatenated.filter { it[1].endsWith(".fastq") } 
    path_ch = concatenated.filter { it[1].contains("/indexed/") }
    // Group fastq and path files by key
    grouped_fastq_ch = fastq_ch.groupTuple(by: 0).map { key, files -> [key, files.flatten()] }
    grouped_path_ch = path_ch.groupTuple(by: 0).map { key, paths -> [key, paths.flatten()] }
    // Combine these channels to get all combinations of fastq and path for each key
    left_joined_ch = grouped_fastq_ch
          .combine(grouped_path_ch, by: 0)
          .flatMap { key, fastqFiles, pathFiles -> 
                 fastqFiles.collect { fastq -> tuple(fastq, pathFiles[0]) }
          }
    typing_input = left_joined_ch.combine(demultiplexed).map { primer, indexed, demultiplexed_dir -> 
        return [primer, indexed, demultiplexed_dir, params.outdir]
    } 
    // | view()

    // TYPING
    typing_res = Typing(typing_input)

    // EXTRACTING 3 BEST MATCHES FOR EACH LOCUS AND COMBINING WITH THE FULL INDEXED DATABASE 
    mapped = typing_res.map { typed -> tuple(typed[0], typed[1]) }
    sorted_typing_res = SortTypingResults(mapped)
    typing_res_from_comparison = CompareProportions(sorted_typing_res)
    selected_records = ExtractNames(typing_res_from_comparison)
    .map{ file -> tuple(getGeneTypeFromName(file.baseName), file) }

    extraction_input = selected_records
    .join(typing_res)
    .map{it.getAt([0,1,3])}

    indexed_hla_path = indexed_ch
    .filter { it[0] == 'hla' } // Filters to only include items where geneType is 'hla'
    .map { it[1] } // Maps the channel to get only the indexPath

    matching_sequences = ExtractMatchingSequences(extraction_input) 
    .map{ file -> tuple(getGeneTypeFromName(file.baseName), file) }
    .combine(indexed_hla_path) // Combining with the HLA index path

    // MAPPING THE BEST RESULTS FROM KMA AGAINST THE WHOLE DATABASE
    final_mapping_res = FinalMapping(matching_sequences)
    results_ch = GatherResults(final_mapping_res.collect())

    // CONVERTING THE RESULTS TO A CSV FILE
    ConvertToCsv(results_ch) | view()
}