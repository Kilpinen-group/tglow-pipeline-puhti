#!/usr/bin/env nextflow

process finalize {
    label params.cpr_label
    
    conda params.cpr_conda_env
    container params.tg_container
    
    storeDir "$params.rn_publish_dir/processed_images"
    scratch params.rn_scratch

    input:
        tuple val(key), val(well), val(manifest), path(cell_masks, stageAs: "cell_masks/*"), path(nucl_masks, stageAs: "nucl_masks/*"), val(decon_plates), val(merge_plates), path(registration, stageAs:"registration_tmp/*")
        path images, stageAs: "input_images"
        tuple path(basicpy_files), val(basicpy_string)
        path scaling_file
        path slope_file
        path bias_file
    output:
        tuple val(well), path("${well.relpath}/*.tiff"), emit: processed_output
        path "${well.plate}/channel_indices.tsv", emit: channel_indices
    script:
    
        // Stage the masks so cellprofiler can access them
        cmd =
        """
        mkdir -p ./masks/${well.relpath}
        ln -s \$(pwd)/cell_masks/*  ./masks/${well.relpath}/
        """
        
        if (!nucl_masks[0].name.startsWith("NO_NUCL_MASK")) {
            cmd += "ln -s \$(pwd)/nucl_masks/* ./masks/${well.relpath}/"
        }
        
        // Stage the registration 
        cmd += 
        """
        mkdir -p registration/${well.plate}/${well.row}/
        ln -s \$(pwd)/registration_tmp/* registration/${well.plate}/${well.row}/
        """
        
        // Outputs the proccessed images
        cmd += 
        """
        # Stage files
        stage_cellprofiler.py \
        --input input_images \
        --output ./ \
        --output_format OME_TIFF \
        --well ${well.well} \
        --plate ${well.plate} \
        """
        
        if (merge_plates) {
            cmd += " --plate_merge " + merge_plates.join(" ")
        }
        
        if (params.rn_manifest_registration) {
            cmd += " --registration_dir ./registration"
        }
                
        if (basicpy_string && basicpy_string != "NO_FLATFIELD") {
            cmd += " --flatfields $basicpy_string"
        }
        
        if ((params.rn_manualscale != null | params.rn_autoscale) & scaling_file.name != "NO_SCALE")  {
            cmd += " --scaling_factors $scaling_file"
        }
        
        if (slope_file.name != "NO_SLOPE")  {
            cmd += " --scaling_slope $slope_file"
        }  
                
        if (bias_file.name != "NO_BIAS")  {
            cmd += " --scaling_bias $bias_file"
        }
        
        if (params.rn_max_project | params.rn_hybrid) {
            cmd += " --max_project --no_zstack"
        }
        
        if (params.rn_hybrid & manifest.mask_channels != "none") {
            cmd += " --mask_dir ./masks"
            cmd += " --mask_pattern *_nucl_mask_*_cp_masks.tiff"
            cmd += " --mask_channels ${manifest.mask_channels.collect{ well.plate + "=" + it }.join(' ')}"
        }
        
        if (params.rn_hybrid) {
            cmd += 
            """
            max_project.py \
            --input ./masks \
            --output ./ \
            --well ${well.well} \
            --plate ${well.plate} \
            --pattern *_cell_mask_*_cp_masks.tiff \
            --suffix _cell_mask_d00_ch0_cp_masks.tiff
            """
            
            if (!nucl_masks[0].fileName.name.startsWith("NO_NUCL_MASK")) {
                cmd += 
                """
                max_project.py \
                --input ./masks \
                --output ./ \
                --well ${well.well} \
                --plate ${well.plate} \
                --pattern *_nucl_mask_*_cp_masks.tiff \
                --suffix _nucl_mask_d00_ch0_cp_masks.tiff
                """
            }

        } else {
            cmd += "\nrsync --copy-links ./masks/${well.relpath}/* ./${well.relpath}/"
        }
        
        cmd
    stub:
        """
        mkdir -p ${well.relpath}
        cd ${well.relpath}
        touch ${well.plate}_${well.well}_ch0.tiff
        cd ../../
        touch channel_indices.tsv
        """ 
}

process cellcrops {
    label params.cpr_label
    
    conda params.tg_conda_env
    container params.tg_container

    storeDir "$params.rn_publish_dir/cellcrops"
    scratch params.rn_scratch

    input:
        tuple val(well), val(registration),  path(images, stageAs: "input_images/")
    output:
        tuple val(well), path("${well.relpath}/*.h5"), emit: h5
        tuple val(well), path("${well.relpath}/*.csv"), emit: index
    script:

    cmd = """
    mkdir -p input/${well.relpath}
    ln -s \$(pwd)/input_images/* input/${well.relpath}
    
    run_cellsampling.py \
    --input input \
    --cell_mask_dir input \
    --nucl_mask_dir input \
    --output ./ \
    --plate ${well.plate} \
    --well ${well.well} \
    --max_per_field $params.rn_max_per_field \
    """   
    if (registration != null) {
        cmd += "--ref_channel ${registration.ref_channel} --qry_channels ${registration.qry_channels.join(" ")}"
    }
    
    cmd
}

process index_cellcrops {
    label params.cpr_label
    
    conda params.tg_conda_env
    container params.tg_container
    
    publishDir "$params.rn_publish_dir/cellcrops", mode: "copy"
    
    input:
        val previous_completed
        path input, stageAs: "input_cellcrops"
    output:
        path "cellcrop_index.csv.gz"
    script:
    """
    # Get the header from the first CSV file
    find input_cellcrops/ -name "*.csv" -type f -print0 | head -z -n1 | xargs -0 head -n1 > cellcrop_index.csv
    
    # Append all CSV files without their headers
    find input_cellcrops/ -name "*.csv" -type f -print0 | xargs -0 tail -q -n+2 >> cellcrop_index.csv
    
    gzip -f cellcrop_index.csv    
    """
}