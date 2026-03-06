#!/usr/bin/env nextflow

include { calculate_scaling_factors; calculate_plate_offsets } from '../processes/scaling.nf'


workflow estimate_scaling_factors {
    
    take:
        manifest
        manifest_registration
        blacklist_file
        control_file
        plates
        manifest_registration_file
        cellpose_out
        flatfield_out
        rn_autoscale
        rn_manualscale
        rn_manifest_registration
        rn_scale_slope
        rn_scale_bias
    main:
    
        // Get scaling string
        if (rn_manualscale != null) {    
            scaling_channel = channel.value(file(rn_manualscale))
        } else {
            scaling_channel = Channel.value(file("NO_SCALE"))
        }
        
        // Optional slope for the sigmoid curve
        if (rn_scale_slope != null) {    
            slope_channel = channel.value(file(rn_scale_slope))
        } else {
            slope_channel = Channel.value(file("NO_SLOPE"))
        }
        
        // Optional bias for the sigmoid curve
        if (rn_scale_bias != null) {    
            bias_channel = channel.value(file(rn_scale_bias))
        } else {
            bias_channel = Channel.value(file("NO_BIAS"))
        }
        
        
        // When all deconvelution is done, or all data is staged, calculate the scaling factors
        // which map the images onto the full dynamic range of the 16 bit uint
        // Alternatively, if control scaling is provided, wait for cellpose and run the controls.
        if (rn_autoscale) {
            // Channel <plate> <mask_channel1 mask_channel2 mask_channelN>
            scaling_in = manifest.map{
                row -> tuple(row.plate,
                (row.mask_channels == "none") ? "none" : row.mask_channels.collect{it -> (row.plate + "=" + it.toString())}.join(" "),
                null)
            }
            
            // Filter to reference plates only as we want scaling factors to be in final dimension
            if (rn_manifest_registration != null){
                scaling_in = scaling_in.combine(manifest_registration, by: 0).map{row -> tuple(row[0], row[1], row[4].split(',').join(" "))}
            }
            
            if (params.rn_control_list) {
                // Start only if registration & cellpose is done
                control_dir = calculate_plate_offsets(
                    cellpose_out.last(), // ensures this only starts when cellpose is done
                    scaling_in,
                    Channel.value(file("${params.rn_publish_dir}/registration")),
                    Channel.value(file("${params.rn_publish_dir}/masks")),
                    flatfield_out,
                    blacklist_file,
                    control_file
                ).plate_offset
            } else {
                control_dir = Channel.value(file('NO_CONTROL_DIR'))
            }
            
            // Final scaling
            demultiplex_channelstring = scaling_in.map{row -> row[1]}.collect().map({it.unique().join(" ")})

            // Start running only if all the plate offsets have been calculated
            scaling_channel = calculate_scaling_factors(cellpose_out.last(),
                blacklist_file,
                plates,
                manifest_registration_file,
                control_dir.collect(), // ensures this only runs when control dir is done
                demultiplex_channelstring).scaling_factors.first() // Emit the scale factor file as a value channel
        
        }
    emit:
        scaling_file=scaling_channel
        slope_file=slope_channel
        bias_file=bias_channel
}