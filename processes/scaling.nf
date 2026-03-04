#!/usr/bin/env nextflow

// Determine scaling factors
process calculate_scaling_factors {
    label 'normal'
    
    conda params.tg_conda_env
    container params.tg_container
    
    storeDir "$params.rn_publish_dir/scaling"
   // publishDir "$params.rn_publish_dir/", mode: 'copy'

    input:
        val x
        path blacklist
        val plates
        path registration_manifest
        path control_dir, stageAs: "control_intensities/*"
        val mask_channels
    output:
        path "scaling_factors.txt", emit: scaling_factors
        path "raw_scaling_factors.txt"
        path "intensity_summary.tsv"
        path "channel_index_with_scaling.tsv"
        path "histograms"
    script:
        cmd = 
        """
        calculate_scaling_factors.py \
        --output ./ \
        --plate $plates \
        --q1 $params.rn_autoscale_q1 \
        --q2 $params.rn_autoscale_q2\
        """
        
        if (params.dc_run) {
            // Not sure why I thought the scale_max was needed. The intensity files match the 16 bit already, so this
            // is not needed.
            //cmd += " --input $params.rn_decon_dir --scale_max $params.dc_clip_max"
            cmd += " --input $params.rn_decon_dir --scale_max 65535"
        } else {
            cmd += " --input $params.rn_image_dir --scale_max 65535"
        }
        
        // Add optional blacklist
        if (params.rn_blacklist) {
            cmd += " --blacklist $blacklist"
        }
        
        if (control_dir.fileName.name != "NO_CONTROL_DIR") {
             cmd += " --control_dir control_intensities"
        }
        
        // add optional registration manifest
        if (params.rn_manifest_registration) {
            cmd += " --plate_groups $registration_manifest"
        }
        
        if (params.rn_hybrid & mask_channels != "none") {
            cmd += " --mask_channels $mask_channels"
        }
        
        // TMP dummy variable
        //cmd = "echo test > scaling_factors.txt"   
        cmd
    stub:
        """
        touch scaling_factors.txt
        touch raw_scaling_factors.txt
        touch intensity_summary.tsv
        touch channel_index_with_scaling.tsv
        mkdir histograms
        """
}

// Determine offsets for scaling factors based on controls
process calculate_plate_offsets {
    label 'normal_plus'
    
    conda params.tg_conda_env
    container params.tg_container
    
    storeDir "$params.rn_publish_dir/scaling/offsets/"
    stageInMode 'symlink'
   // publishDir "$params.rn_publish_dir/", mode: 'copy'
   
    input:
        val x
        tuple val(plate), val(mask_channels), val(merge_plates)
        path registration, stageAs:"registration"
        path masks, stageAs: 'masks'
        tuple path(bp_files), val(basicpy_string)
        path blacklist
        path control_list
    output:
        path "$plate", emit: plate_offset
    script:
        cmd = 
        """
        masked_control_intensity_calculator.py \
        --output ./ \
        --plate $plate \
        --obj_mask_dir ./masks \
        --obj_mask_pattern *_cell_mask_*_cp_masks.tiff\
        """
        
        if (params.rn_dummy_mode) {
            cmd += " --dummy_mode"
        }
        
        if (params.rn_threshold) {
            cmd += " --threshold"
        }
        
        if (params.rn_control_list) {
            cmd += " --controls $control_list"
        }
        
        if (params.rn_blacklist) {
            cmd += " --blacklist $blacklist"
        }
        
        if (params.dc_run) {
            cmd += " --input $params.rn_decon_dir"
        } else {
            cmd += " --input $params.rn_image_dir"
        }
        
        if (merge_plates) {
            cmd += " --plate_merge " + merge_plates
        }
        
        if (registration.fileName.name != "NO_REGISTRATION") {
            cmd += " --registration_dir ./registration"
        }
                
        if (basicpy_string && basicpy_string != "NO_FLATFIELD") {
            cmd += " --flatfields $basicpy_string"
        }
        
        if (params.rn_max_project | params.rn_hybrid) {
            cmd += " --max_project"
        }
        
        if (params.rn_hybrid & mask_channels != "none") {
            cmd += " --mask_dir ./masks"
            cmd += " --mask_pattern *_nucl_mask_*_cp_masks.tiff"
            cmd += " --mask_channels $mask_channels"
        }
    stub:
        """
        mkdir $plate
        cd $plate
        touch plate_offset.txt
        """

}
