#!/usr/bin/env nextflow

params.output = "./output/"
params.meta_file = "./example/metadata.json"
params.count_files = "./example/data/00/00/03/1/quantification_table.tsv,./example/data/00/00/04/1/quantification_table.tsv"

metadata = Channel.fromPath(params.meta_file)
file_list = params.count_files.tokenize(",")
file_channels = Channel.fromPath(file_list).collect()

// scripts
deanalysis_script = Channel.fromPath("${projectDir}/docker_deanalysis/DEAnalysis.R")
join_table = Channel.fromPath("${projectDir}/semares_preprocessing/join_table.py")
metadata2table = Channel.fromPath("${projectDir}/semares_preprocessing/metadata2table.py")

// config files
data_config = Channel.fromPath("${projectDir}/config/data_table_config.json")
meta_data_config = Channel.fromPath("${projectDir}/config/meta_table_config.json")

process file_join {
    container "dockergenevention/pandas" // use docker conatainer
    publishDir params.output, mode: "copy"

    input:
    path script
    path metadata
    path config
    path file_channels, stageAs: "*/*"
    val file_path

    output:
    path "output/data.tsv"

    """
    python $script -m $metadata -c $config -f $file_channels -p $file_path
    """
}

process metadata_join {
    container "dockergenevention/pandas" // use docker conatainer
    publishDir params.output, mode: "copy"

    input:
    path script
    path metadata
    path config
    val file_path

    output:
    path "output/metadata.tsv"

    """
    python $script -m $metadata -c $config -p $file_path
    """
}



process deanalysis {
    container 'kadam0/deanalysis:0.0.1'
    publishDir params.output, mode: "copy"

    input:
    path script_file
    path meta_file 
    path count_file 

    output:                                
    path "*"

    script:
    """
    Rscript $script_file --meta_file ${meta_file} --count_file ${count_file} --out_dir ./
    """
}

workflow {
  file_join(join_table, metadata, data_config, file_channels, params.count_files)
  metadata_join(metadata2table, metadata, meta_data_config, params.count_files)
  deanalysis(deanalysis_script, metadata_join.out, file_join.out)
}
