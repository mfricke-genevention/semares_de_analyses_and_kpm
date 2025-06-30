#!/usr/bin/env nextflow

params.output = "./output/"
params.meta_file = "./example/metadata.json"
params.count_files = "./example/data/1/quantification_table.tsv,./example/data/2/quantification_table.tsv,./example/data/3/quantification_table.tsv,./example/data/4/quantification_table.tsv,./example/data/5/quantification_table.tsv,./example/data/6/quantification_table.tsv,./example/data/7/quantification_table.tsv,./example/data/8/quantification_table.tsv,./example/data/10/quantification_table.tsv"

params.logFC = true
params.logFC_up = 1
params.logFC_down = -1
params.p_adj = true
params.alpha = 0.05
params.gene_column = "RemappedGeneNames"
params.organism = "rat"

params.de_analysis = false
params.kpm_analysis = false
params.go_enrichment = false
params.network_enrichment = true


params.file_type = "png"
params.protein_column = "Protein.IDs"
params.condition_1 = ""
params.condition_2 = ""
params.condition_3 = ""
params.rev_con = false
params.reviewed = false
params.mode = "uniprot_one"
params.skip_filled = true
params.fasta = null
params.tar_organism = "human"


metadata = Channel.fromPath(params.meta_file)
file_list = params.count_files.tokenize(",")
file_channels = Channel.fromPath(file_list).collect()

network_file = Channel.fromPath("${projectDir}/network.sif")

// scripts
deanalysis_script = Channel.fromPath("${projectDir}/docker_deanalysis/DEAnalysis.R")
summarizer_script = Channel.fromPath("${projectDir}/docker_deanalysis/Summarizer.R")
kpm_script = Channel.fromPath("${projectDir}/docker_kpm/KPMAnalysis.R")
go_enrichment_script = Channel.fromPath("${projectDir}/docker_goenrichment/GOEnrichment.R")
network_enrichment_script = Channel.fromPath("${projectDir}/docker_networkenrichment/NetworkEnrichment.R")
join_table = Channel.fromPath("${projectDir}/semares_preprocessing/join_table.py")
metadata2table = Channel.fromPath("${projectDir}/semares_preprocessing/metadata2table.py")
update_metadata_config = Channel.fromPath("${projectDir}/update_metadata_config.py")
proharmed_script = Channel.fromPath("${projectDir}/docker_proharmed/ProHarMeD.R")
proharmed_config_script = Channel.fromPath("${projectDir}/docker_proharmed/config.R")


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
    path update_metadata_config_script
    path metadata
    path config
    val file_path

    output:
    path "transformed_input/metadata.tsv"

    """
    python $update_metadata_config_script -a $params.condition_1 -b $params.condition_2 -c $params.condition_3 -m $config
    python $script -m $metadata -c ./meta_table_config_changed.json  -p $file_path
    """
}


process proharmed {
    container "kadam0/proharmed:0.0.3"
    publishDir params.output, mode: "copy"

    input:
    path script_file
    path config_file
    path count_file 
    val file_type
    val protein_column
    val organism
    val rev_con
    val reviewed
    val mode
    val skip_filled
    val tar_organism

    output:                                
    path("proharmed_output", type:"dir")
    path("proharmed_output/harmonized_data.csv", emit: "harmonized_data")

    script:
    """
    Rscript ${script_file} --count_file ${count_file} --out_dir ./proharmed_output --file_type ${file_type} --protein_column ${protein_column} --organism ${organism} --rev_con ${rev_con} --reviewed ${reviewed} --mode ${mode}  --skip_filled ${skip_filled} --tar_organism ${tar_organism}
    """
}


process deanalysis {
    container "kadam0/deanalysis:0.0.2"
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
    container "kadam0/kpmanalysis:0.0.2"
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
    val gene_column
    val organism

    output:
    path("go_enrichment_out", type:"dir")

    """
    Rscript $go_enrichment_script --meta_file ${meta_file} --count_file ${count_file} --out_dir ./go_enrichment_out --logFC ${logFC} --logFC_up ${logFC_up} --logFC_down ${logFC_down} --p_adj ${p_adj} --alpha ${alpha} --gene_column ${gene_column} --organism ${organism}
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
    val gene_column

    output:
    path("network_enrichment_out", type:"dir")

    """
    Rscript $network_enrichment_script --meta_file ${meta_file} --count_file ${count_file} --out_dir ./network_enrichment_out --logFC ${logFC} --logFC_up ${logFC_up} --logFC_down ${logFC_down} --p_adj ${p_adj} --alpha ${alpha} --gene_column ${gene_column}
    """
}

workflow {
  file_join(join_table, metadata, data_config, file_channels, params.count_files)
  metadata_join(metadata2table, update_metadata_config, metadata, meta_data_config, params.count_files)
//   proharmed(proharmed_script, proharmed_config_script, file_join.out, params.file_type, params.protein_column, params.organism, params.rev_con, params.reviewed, params.mode, params.skip_filled, params.tar_organism)

//   if( params.de_analysis ){
//     deanalysis(deanalysis_script, metadata_join.out, file_join.out)
//     summarize(summarizer_script, deanalysis.out, params.logFC, params.logFC_up, params.logFC_down, params.p_adj, params.alpha)
//   }
//   if( params.kpm_analysis )
//     kpm_analysis(kpm_script, file_join.out, metadata_join.out, network_file, params.logFC, params.logFC_up, params.logFC_down, params.p_adj, params.alpha)
//   if( params.go_enrichment )
//     go_enrichment(go_enrichment_script, proharmed.out.harmonized_data, metadata_join.out, network_file, params.logFC, params.logFC_up, params.logFC_down, params.p_adj, params.alpha, params.gene_column, params.organism)
//   if( params.network_enrichment )
//     network_enrichment(network_enrichment_script, proharmed.out.harmonized_data, metadata_join.out, network_file, params.logFC, params.logFC_up, params.logFC_down, params.p_adj, params.alpha, params.gene_column)
}
