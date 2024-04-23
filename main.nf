#!/usr/bin/env nextflow

params.output = "./output/"
params.meta_file = "./example/metadata.json"
params.count_files = "./example/data/1/quantification_table.tsv,./example/data/2/quantification_table.tsv,./example/data/3/quantification_table.tsv,./example/data/4/quantification_table.tsv,./example/data/5/quantification_table.tsv,./example/data/6/quantification_table.tsv,./example/data/7/quantification_table.tsv,./example/data/8/quantification_table.tsv,./example/data/10/quantification_table.tsv"

params.logFC = true
params.logFC_up = 1
params.logFC_down = -1
params.p_adj = true
params.alpha = 0.05


metadata = Channel.fromPath(params.meta_file)
file_list = params.count_files.tokenize(",")
file_channels = Channel.fromPath(file_list).collect()

network_file = Channel.fromPath("${projectDir}/network.sif")

// scripts
deanalysis_script = Channel.fromPath("${projectDir}/docker_deanalysis/DEAnalysis.R")
summarizer_script = Channel.fromPath("${projectDir}/docker_deanalysis/Summarizer.R")
kpm_script = Channel.fromPath("${projectDir}/docker_kpm/KPMAnalysis.R")
go_enrichment_script = Channel.fromPath("${projectDir}/docker_geoenrichment/GOEnrichment.R")
network_enrichment_script = Channel.fromPath("${projectDir}/docker_networkenrichment/NetworkEnrichment.R")
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
    path "transformed_input/data.tsv"

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
    path "transformed_input/metadata.tsv"

    """
    python $script -m $metadata -c $config -p $file_path
    """
}



process deanalysis {
    container 'kadam0/deanalysis:0.0.2'
    publishDir params.output, mode: "copy"

    input:
    path script_file
    path meta_file 
    path count_file 

    output:                                
    path("de_out", type:"dir")

    script:
    """
    Rscript $script_file --meta_file ${meta_file} --count_file ${count_file} --out_dir ./de_out
    """
}


process summarize {
    container 'kadam0/kpmanalysis:0.0.2'
    publishDir params.output, mode: "copy"

    input:
    path script_file
    path input_dir
    val logFC
    val logFC_up
    val logFC_down
    val p_adj
    val alpha

    output:                                
    path "summary"

    script:
    """
    Rscript $script_file --in_dir ${input_dir} --out_dir ./summary --logFC ${logFC} --logFC_up ${logFC_up} --logFC_down ${logFC_down} --p_adj ${p_adj} --alpha ${alpha}    """
}

process kpm_analysis {
    container "kadam0/kpmanalysis:0.0.2" // use docker conatainer
    publishDir params.output, mode: "copy"

    input:
    path kpm_script
    path count_file
    path meta_file
    path network_file
    val logFC
    val logFC_up
    val logFC_down
    val p_adj
    val alpha

    output:
    path("kpm_out", type:"dir")

    """
    Rscript $kpm_script --meta_file ${meta_file} --count_file ${count_file} --network_file ${network_file} --out_dir ./kpm_out --logFC ${logFC} --logFC_up ${logFC_up} --logFC_down ${logFC_down} --p_adj ${p_adj} --alpha ${alpha}
    """
}

process go_enrichment {
    container "kadam0/goenrichment:0.0.2" // use docker conatainer
    publishDir params.output, mode: "copy"

    input:
    path go_enrichment_script
    path count_file
    path meta_file
    path network_file
    val logFC
    val logFC_up
    val logFC_down
    val p_adj
    val alpha

    output:
    path("go_enrichment_out", type:"dir")

    """
    Rscript $go_enrichment_script --meta_file ${meta_file} --count_file ${count_file} --network_file ${network_file} --out_dir ./go_enrichment_out --logFC ${logFC} --logFC_up ${logFC_up} --logFC_down ${logFC_down} --p_adj ${p_adj} --alpha ${alpha}
    """
}

process  network_enrichment {
    container "kadam0/netenrichment:0.0.2" // use docker conatainer
    publishDir params.output, mode: "copy"

    input:
    path network_enrichment_script
    path count_file
    path meta_file
    path network_file
    val logFC
    val logFC_up
    val logFC_down
    val p_adj
    val alpha

    output:
    path("network_enrichment_out", type:"dir")

    """
    Rscript $network_enrichment_script --meta_file ${meta_file} --count_file ${count_file} --network_file ${network_file} --out_dir ./network_enrichment_out --logFC ${logFC} --logFC_up ${logFC_up} --logFC_down ${logFC_down} --p_adj ${p_adj} --alpha ${alpha}
    """
}

workflow {
  file_join(join_table, metadata, data_config, file_channels, params.count_files)
  metadata_join(metadata2table, metadata, meta_data_config, params.count_files)
  deanalysis(deanalysis_script, metadata_join.out, file_join.out)
  summarize(summarizer_script, deanalysis.out, params.logFC, params.logFC_up, params.logFC_down, params.p_adj, params.alpha)
  kpm_analysis(kpm_script, file_join.out, metadata_join.out, network_file, params.logFC, params.logFC_up, params.logFC_down, params.p_adj, params.alpha)
  go_enrichment(kpm_script, file_join.out, metadata_join.out, network_file, params.logFC, params.logFC_up, params.logFC_down, params.p_adj, params.alpha)
  network_enrichment(kpm_script, file_join.out, metadata_join.out, network_file, params.logFC, params.logFC_up, params.logFC_down, params.p_adj, params.alpha)
}
