#!/usr/bin/env nextflow


// Cellpose
process cellpose {
    label params.cp_label
    
    conda params.tg_conda_env
    container params.tg_container
    
    storeDir "${params.rn_publish_dir}/masks/"
    scratch params.rn_scratch

    input:
        tuple val(well), val(manifest), path(images, stageAs: "input_images/")
    output:
        tuple val(well), val(manifest), path("${well.relpath}/*_cell_mask*_ch${manifest.cp_cell_channel}*.*"), path("${well.relpath}/*_nucl_mask*_ch${manifest.cp_nucl_channel}*.*"), emit: cellpose_out
    script:
        cmd =
        """
        # Workarround as we cannot use variables from the same tuple in stageAs
        mkdir -p input/${well.plate}/${well.row}
        ln -s \$(pwd)/input_images/* input/${well.plate}/${well.row}
        
        run_cellpose.py \
        --input input/ \
        --output ./ \
        --plate ${well.plate} \
        --well ${well.well} \
        --cell_channel ${manifest.cp_cell_channel} \
        --gpu \
        --diameter ${params.cp_cell_size} \
        --model ${params.cp_model}\
        """

        if (manifest.cp_nucl_channel >= 0) {
            cmd +=
            """ \
            --nucl_channel ${manifest.cp_nucl_channel} \
            --diameter_nucl ${params.cp_nucl_size} \
            """
        }    
        
        if (params.cp_min_cell_area) {
            cmd += " --min_cell_area ${params.cp_min_cell_area}"
        }
        
        if (params.cp_min_nucl_area) {
            cmd += " --min_nucl_area ${params.cp_min_nucl_area}"
        }
        
        if (params.cp_cell_flow_thresh) {
            cmd += " --cell_flow_thresh ${params.cp_cell_flow_thresh}"
        }
        
        if (params.cp_nucl_flow_thresh) {
            cmd += " --nucl_flow_thresh ${params.cp_nucl_flow_thresh}"
        }
        
        if (params.cp_cell_prob_threshold) {
            cmd += " --cell_prob_threshold ${params.cp_cell_prob_threshold}"
        }
        
        if (params.cp_nucl_prob_threshold) {
            cmd += " --nucl_prob_threshold ${params.cp_nucl_prob_threshold}"
        }
        
        if (params.rn_max_project || params.rn_hybrid) {
            cmd += " --no_3d"
        }
        
        if (params.cp_dont_use_nucl_for_declump) {
            cmd += " --dont_use_nucl_for_declump"
        }
        
        if (params.cp_cell_power) {
            cmd += " --cell_power ${params.cp_cell_power}"
        }
        
        if (params.cp_nucl_power) {
            cmd += " --nucl_power ${params.cp_nucl_power}"
        }
        
        if (params.cp_dont_post_process) {
            cmd += " --dont_post_process"
        }
        
        if (params.cp_downsample) {
            cmd += " --downsample ${params.cp_downsample}"
        }
        
        // Add a fake nucleus channel because nextflow doesn't play nicely with 
        // tuples and mulptiple input files, otherwise making sure wells and channels are 
        // matched turns into a pain
        // If its stupid and it works, it is not stupid, hopefully in future, nextflow
        // will deal better with optional files or just accept null for file objects
        if (manifest.cp_nucl_channel < 0) {
            cmd += 
            """
            
            touch ${well.relpath}/NO_NUCL_MASK_nucl_mask_ch${manifest.cp_nucl_channel}_dummy.tif
            """
        
        }
        cmd
    stub:
        """
        mkdir -p "${well.relpath}"
        cd "${well.relpath}"
        touch 1_cell_mask_d38_ch${manifest.cp_cell_channel}_cp_masks.tiff
        touch 1_nucl_mask_d30_ch${manifest.cp_nucl_channel}_cp_masks.tiff
        """
}
