#!/usr/bin/env nextflow

// Run a cellprofiler run
// regular queue
process finalize_and_cellprofiler {
    label params.cpr_label
    
    conda params.cpr_conda_env
    container params.cpr_container
    
    storeDir "$params.rn_publish_dir/features/cellprofiler"
    scratch params.rn_scratch

    input:
        tuple val(key), val(well), val(manifest), path(cell_masks, stageAs: "cell_masks/*"), path(nucl_masks, stageAs: "nucl_masks/*"), val(decon_plates), val(merge_plates), path(registration, stageAs:"registration_tmp/*")
        path images, stageAs: "input_images"
        tuple path(basicpy_files), val(basicpy_string)
        path scaling_file
        path slope_file
        path bias_file
        path pipeline

    output:
        path "features/${well.relpath}/*"
    script:
    
        // Stage the masks
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
    
        // Outputs the cp files into ./images
        cmd += 
        """
        # Stage files
        stage_cellprofiler.py \
        --input ./input_images \
        --output ./images \
        --output_format CP \
        --well ${well.well} \
        --plate ${well.plate}\
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
            --output ./images \
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
                --output ./images \
                --well ${well.well} \
                --plate ${well.plate} \
                --pattern *_nucl_mask_*_cp_masks.tiff \
                --suffix _nucl_mask_d00_ch0_cp_masks.tiff
                """
            }

        } else {
            cmd += "\nln -s ./masks/${well.relpath}/* ./images/${well.relpath}/"
        }
    
        // Run cell profiler
        cmd +=
        """
        # Run cellprofiler
        cellprofiler \
        -p ${pipeline} \
        -c \
        -r \
        -o ./features/${well.relpath} \
        -i ./images/${well.relpath}\
        """
        
        if (params.cpr_plugins) {
            cmd += " --plugins-directory $params.cpr_plugins"
        }
        
        // Zip output to save of lustre filecount
        if (!params.cpr_no_zip) {
            cmd +=
            """
            zip -r ./features/${well.relpath}/${well.plate}_${well.well}.zip ./features/${well.relpath}/*.txt

            # Cleanup so only zip is staged
            rm ./features/${well.relpath}/*.txt
            """
        }

        cmd
    stub:
        """
        mkdir -p features/${well.relpath}
        cd features/${well.relpath}
        touch ${well.plate}_${well.well}.txt
        """
}

process cellprofiler {
    label params.cpr_label
    
    conda params.cpr_conda_env
    container params.cpr_container
    
    storeDir "$params.rn_publish_dir/features/cellprofiler"
    scratch params.rn_scratch

    input:
        tuple val(well), path(images, stageAs: "input_images/*")
        path pipeline
    output:
        path "features/${well.relpath}/*"
    script:
    
        // Outputs the cp files into ./images
        cmd = 
        """
        mkdir -p images/${well.relpath}
        ln -s \$(pwd)/input_images/* images/${well.relpath}/
        
        # Convert format from OME tiff to CP
        stage_cellprofiler.py \
        --input ./images \
        --output ./images \
        --output_format CP \
        --well ${well.well} \
        --plate ${well.plate}
        """
        
        // Run cell profiler
        cmd +=
        """
        # Run cellprofiler
        cellprofiler \
        -p ${pipeline} \
        -c \
        -r \
        -o ./features/${well.relpath} \
        -i ./images/${well.relpath} \
        """
        
        if (params.cpr_plugins) {
            cmd += " --plugins-directory $params.cpr_plugins"
        }
    
        // Zip output to save of lustre filecount
        if (!params.cpr_no_zip) {
            cmd +=
            """
            zip -r ./features/${well.relpath}/${well.plate}_${well.well}.zip ./features/${well.relpath}/*.txt

            # Cleanup so only zip is staged
            rm ./features/${well.relpath}/*.txt
            """
        }

        cmd
    stub:
        """
        mkdir -p features/${well.relpath}
        cd features/${well.relpath}
        touch ${well.plate}_${well.well}.txt
        """ 
    
}