#' @name create_somatic_json
#' @title create_somatic_json
#' @description
#'
#' @param somatic_snv_cn somatic copy number data table
#' @param out_file where the json should be written
#' @param pair patient.id to be used for pgvdb or case reports
#' @param pgv_settings the settings file used for pgvdb or case reports
#' @param return_table_pgv TRUE/FALSE whether to return the data.table to add to pgvdb
#' @param meta_keep meta data to write in the json
#' @param y_col column to be used for y axis values
#' @return data.table to add to pgv or NULL, depends on return_table_pgv
#' @export
#' @author Stanley Clarke
create_somatic_json = function(somatic_snv_cn, out_file, pair, pgv_settings, return_table_pgv = FALSE, meta_keep = NULL, y_col = "est_cn_llrm") {
    som.dt = readRDS(somatic_snv_cn)
    som.dt = som.dt[!is.na(get(y_col)),]
    som.dt[start == end, end := end +1]
    som.dt[, strand := NULL]
    som.dt[variant.p != "",annotation := paste0("Type: ", annotation, "; Gene: ", gene, "; Variant: ",variant.c, "; Protein_variant: ", variant.p, "; VAF: ",vaf)]
    som.dt[variant.p == "",annotation := paste0("Type: ", annotation, "; Gene: ", gene, "; Variant: ",variant.c, "; VAF: ",vaf)]
    dt2json_mut(dt = som.dt, ref = "hg19",settings = pgv_settings, meta_data = meta_keep, y_col = y_col,
                file_name = out_file)
    if(return_table_pgv) {
        dt.add = data.table(patient.id = pair, type = "genome",visible = TRUE, title = "Copy Number Mutations", source = "mutations.json")
    }
}

#' @name dt2json_mut
#' @title dt2json_mut
#' @description
#'
#' function to create a mutation json, used in create_somatic_json
#' 
#' @param dt data.table with seqnames,start and end
#' @param patient.id patient id to be added to pgvdb or case reports
#' @param ref reference for pgv or case reports
#' @param file_name the file the json should be written to
#' @param meta_data extra information to add to the json
#' @param y_col column to be used for y axis values
#' @return NULL
#' @export
#' @author Stanley Clarke
dt2json_mut = function(dt,patient.id,ref,settings,file_name = paste(getwd(),"test.json",sep="/"), meta_data = NULL, y_col = NULL) {
    #create vector of seqlengths
    settings_data <- jsonlite::fromJSON(settings)
    chrom_lengths <- as.data.table(settings_data$coordinates$sets[[ref]])[,.(chromosome,startPoint,endPoint)]
    colnames(chrom_lengths) = c("seqnames","start","end")

    if(nrow(chrom_lengths[grepl("chr",seqnames),]) > 0) {
        chrom_lengths[!grepl("chr",seqnames), seqnames := paste0("chr",seqnames)] # weird fix because hg38_chr does not have chr on Y and M
    }
                                        #add y value specified
    if(is.null(y_col)) {
        dt$y_value = 1
    } else {
        dt[,y_value := get(y_col)]
    }
#convert to ggraph and create json
    if(nrow(chrom_lengths[grepl("chr",seqnames),]) > 0) {
        gr1 = dt2gr(dt[order(seqnames,start),]) %>% sortSeqlevels() %>% gr.chr()
    } else {
        gr1 = dt2gr(dt[order(seqnames,start),]) %>% sortSeqlevels() %>% gr.nochr()
    }
    if(any(gr1@seqinfo@seqlengths > chrom_lengths[seqnames %in% names(seqlengths(gr1))]$end)) {
        stop(paste("the seqlengths of your granges has ranges that are not contained in the seqlengths of",ref))
    }
    jab = gG(nodes=gr1)
    settings_y = list(y_axis = list(title = "copy number",
                                    visible = TRUE))
    node.json = gr2dt(jab$nodes$gr[, c("snode.id","y_value",meta_data)])[, .(chromosome = seqnames, startPoint = start, endPoint = end, iid = snode.id,title=snode.id,type="interval", y = y_value, annotation = annotation)]
    gg.js = list(intervals = node.json, connections = data.table())
    gg.js = c(list(settings = settings_y), gg.js)
    message(paste0("Writing json to ",file_name))
    jsonlite::write_json(gg.js, file_name,
                         pretty=TRUE, auto_unbox=TRUE, digits=4)
    return(NULL)
}

#' @name filtered_events_json
#' @title filtered_events_json
#' @description
#'
#' function to create filtered events json for case reports
#' 
#' @param pair patient id to be added to pgvdb or case reports
#' @param oncotable oncotable task output
#' @param jabba_gg JaBbA output ggraph or complex
#' @param out_file path to write json
#' @param cgc_file path to cgc file to annotate drivers
#' @param return_table TRUE/FALSE whether to return the data.table that is used for creating the json
#' @return data.table or NULL
#' @export
#' @author Stanley Clarke, Tanubrata Dey

filtered_events_json = function(pair, oncotable, jabba_gg, out_file, cgc_file = "/gpfs/commons/groups/imielinski_lab/DB/COSMIC/v99_GRCh37/cancer_gene_census_fixed.csv", return_table = FALSE) {
    ##Driver CNA windows
    ##Load details from oncotable
    ot = readRDS(oncotable)
    snvs = ot[grepl('frameshift|missense|stop|disruptive', annotation)]
    snvs = snvs[!duplicated(variant.p)]
    ##Note here probably have to crossreference these missense muts with hetdels
    #hetdel_snvs = snvs[gene %in% ot[type == 'hetdel',gene]]
                                        #possible_drivers = rbind(hetdel_snvs,homdels)
    homdels = ot[type == 'homdel']
    amps = ot[type == 'amp']
    jab = readRDS(jabba_gg)
    possible_drivers = rbind(snvs,homdels,amps)
    cgc = fread(cgc_file)
    names(cgc) = gsub(' ','.', names(cgc))
    cgc$gene = cgc$Gene.Symbol
    ## longlist = merge.data.table(possible_drivers, cgc, by = 'gene', all.x = TRUE)
    longlist = merge.data.table(possible_drivers, cgc, by = 'gene')
    res = longlist[ ,.(gene, id, type, variant.p, Name, Genome.Location, Tier, Role.in.Cancer)]
    names(res) = c("gene", "id", "type", "Variant", "Name", "Genome_Location", "Tier", "Role_in_Cancer")
                                        #add copy number to homdels
    res = res %>% unique
    if(nrow(res) > 0) {
        res[,seqnames := tstrsplit(Genome_Location,":",fixed=TRUE,keep=1)]
        res[,start := tstrsplit(Genome_Location,"-",fixed=TRUE,keep=1)]
        res[,start := tstrsplit(start,":",fixed=TRUE,keep=2)]
        res[,end := tstrsplit(Genome_Location,"-",fixed=TRUE,keep=2)]
        res.mut = res[!is.na(Variant),]
        if(nrow(res.mut) > 0) {
            res.mut[,Variant := gsub("p.","",Variant)]
        }
        res.cn = res[is.na(Variant),]
        if(nrow(res.cn) >0) {
            res.cn.gr = GRanges(res.cn)
            res.cn.gr = gr.val(res.cn.gr,jab$nodes$gr,c("cn","cn.low","cn.high"))
            res.cn.dt = as.data.table(res.cn.gr)
            res.cn.dt[!is.na(cn) & !is.na(cn.low) & !is.na(cn.high), Variant := paste0("Total CN:",round(cn,digits = 3),"; CN Low:",round(cn.low,digits = 3),"; CN High:",round(cn.high,digits = 3))]
            res.cn.dt[!is.na(cn) & is.na(cn.low) & is.na(cn.high), Variant := paste0("Total CN:",round(cn,digits = 3))]
            res.cn.dt[,cn := NULL]
            res.cn.dt[,cn.high := NULL]
            res.cn.dt[,cn.low := NULL]
            res.cn.dt[,width := NULL]
            res.cn.dt[,strand := NULL]
            res.final = rbind(res.mut,res.cn.dt)
        } else {
            res.final = res.mut
            res.final[,seqnames := NULL]
            res.final[,start := NULL]
            res.final[,end := NULL]
        }
        write_json(res.final, out_file, pretty=TRUE)
        res.final[,sample := pair]
        if(return_table) {
            return(res.final)
        } else {
            return(NULL)
        }
    }
}

#' @name dlrs
#' @title derivative log ratio spread
#' @description
#'
#' function to get the dlrs for coverage data
#' used in meta_data_json
#' 
#' @param x foreground from dryclean
#' @return dlrs
#' @export
#' @author Joel Rosiene
dlrs = function(x) {
    nx = length(x)
    if (nx<3) {
        stop("Vector length>2 needed for computation")
    }
    tmp = embed(x,2)
    diffs = tmp[,2]-tmp[,1]
    dlrs = IQR(diffs,na.rm=TRUE)/(sqrt(2)*1.34)
    return(dlrs)
}


#' @name meta_data_json
#' @title meta_data_json
#' @description
#'
#' function to create the meta data summary json for case reports
#' 
#' @param pair patient id to be added to pgvdb or case reports
#' @param out_file path to write json
#' @param coverage path dryclean coverage output
#' @param jabba_gg path to JaBbA output ggraph or complex
#' @param vcf path to strelka vcf to get snv count
#' @param svaba_somatic_vcf path to svaba somatic vcf for getting sv count
#' @param tumor_type_final tumor type abbreviation of the sample
#' @param disease full length tumor type
#' @param primary_site primary site of tumor
#' @param inferred_sex sex of the patient
#' @param karyograph JaBbA outputted karygraph
#' @param seqnames_loh chromosomes to be used to calculate LOH
#' @param seqnames_genome_width chromosomes to be used to calculate tmb
#' @param write_json TRUE/FALSE to write the json
#' @param overwrite TRUE/FALSE to overwrite the present json
#' @param return_table TRUE/FALSE to return the data.table output
#' @return data.table or NULL
#' @export
#' @author Stanley Clarke, Tanubrata Dey, Joel Rosiene

meta_data_json = function(pair, out_file, coverage, jabba_gg, vcf, svaba_somatic_vcf, tumor_type_final, disease, primary_site, inferred_sex, karyograph, seqnames_loh = c(1:22), seqnames_genome_width = c(1:22,"X","Y"), write_json = TRUE, overwrite = FALSE, return_table = FALSE) {
    if(!overwrite) {
        if(file.exists(out_file)) {
            print('Output already exists! - skipping sample')
            return(NA)
        }
    }
    meta.dt = data.table(pair = pair, tumor_type_final = tumor_type_final, disease = disease, primary_site = primary_site, inferred_sex = inferred_sex)
    #get derivate log ratio spread
    meta.dt$dlrs = dlrs(readRDS(coverage)$foreground)
    ##Load this case's counts
    ## meta.dt$snv_count = length(read_vcf(vcf))
    vcf.gr = read_vcf(vcf)
    meta.dt$snv_count = length(vcf.gr %Q% (seqnames %in% seqnames_genome_width))
    ##Count variants
    cmd = paste0("module unload java && module load java; module load gatk; gatk CountVariants --QUIET true --verbosity ERROR"," -V ",svaba_somatic_vcf)
    meta.dt$sv_count = system(paste(cmd, "2>/dev/null"), intern = TRUE)[2] %>% as.integer() #run the command without printing the java command
                                        #Load jabba_gg and get loh
    gg = readRDS(jabba_gg)
    nodes.dt = gg$nodes$dt[seqnames %in% seqnames_loh]
    totalseglen = nodes.dt$width %>% sum()
    if('cn.low' %in% names(nodes.dt)) {
        LOHsegs = nodes.dt[cn.low==0,] %>% .[cn.high >0] %>% .$width %>% sum()
        ## LOH
        LOH_frc = LOHsegs/totalseglen
        meta.dt[,loh_fraction := LOH_frc]
        meta.dt[,loh_seglen := LOHsegs]
        meta.dt[,loh_total_genome := totalseglen]
    } else {
        meta.dt$loh_fraction = 'Not Allelic Jabba'
    }
    #add purity and ploidy
    meta.dt$purity = gg$meta$purity
    meta.dt$ploidy = gg$meta$ploidy

    #' Load beta/gamma for karyograph
    kag = readRDS(karyograph)
    meta.dt$beta = kag$beta
    meta.dt$gamma = kag$gamma
    #add the total seqlengths by using the seqlengths in the jabba object
    nodes.gr = gg$nodes$gr
    seqlengths.dt = seqinfo(nodes.gr) %>% as.data.table(.,keep.rownames = "seqnames")
    seqlengths.dt = seqlengths.dt[seqnames %in% seqnames_genome_width,]
    meta.dt$total_genome_length = sum(seqlengths.dt$seqlengths)
                                        #add tmb
    meta.dt[,tmb := (snv_count / (as.numeric(meta.dt$total_genome_length) / 1e6))]
    if(write_json) {
                                        #write the json
        message(paste0("Writing json to ",out_file))
        write_json(meta.dt,out_file,pretty = TRUE)
        if(return_table) {
            return(meta.dt)
        }
        
    } else {
        return(meta.dt)
    }
}


#' @name strelka_qc
#' @title strelka_qc
#' @description
#'
#' function to create json for strelka qc plotting in case reports
#' 
#' @param strelkaqc_filtered_rds path to an rds of filtered strelka output
#' @param outfile path to write json
#' @param write_json TRUE/FALSE whether to write the json
#' @param return_table TRUE/FALSE whether to return the data
#' @return NULL
#' @export 
#' @author Tanubrata Dey
strelka_qc = function(strelkaqc_filtered_rds, outfile, write_json = TRUE, return_table = TRUE) {
    strelka.qc = readRDS(strelkaqc_filtered_rds) %>% as.data.table
    sq = strelka.qc[ ,.(CHROM,POS,REF,ALT,T_DP,N_DP,MQ,VAF_T, somatic_EVS)]
    names(sq) = c("chromosome", "position", "reference", "alternate", "tumor_depth", "normal_depth", "mapping_quality", "tumor_VAF", "somatic_EVS")
    setnames(sq, old = names(sq), new = new_names_sq)
    if(write_json) {
                                        #write the json
        message(paste0("Writing json to ",outfile))
        write_json(sq,outfile,pretty = TRUE)
        if(return_table) {
            return(sq)
        }

    } else {
        return(sq)
    }
}


#' @name create_distributions
#' @title create_distributions
#' @description
#'
#' function to create the meta data summary for samples in a case reports instance
#' 
#' @param case_reports_data_folder folder with all case report data by sample
#' @param common_folder path to a folder to write all 7 jsons
#' @param filter_pateints list of samples to filter on to create the distributions
#' @return NULL
#' @export
#' @author Stanley Clarke, Tanubrata Dey, Joel Rosiene

create_distributions = function(case_reports_data_folder,common_folder, filter_patients = NULL) {
    files.lst = list.files(case_reports_data_folder)
    files.lst = grep("data",files.lst,invert=TRUE, value = TRUE)
    meta.dt = data.table(meta_json = paste0(case_reports_data_folder,files.lst,"/metadata.json"), patient_id = file.lst)
    meta.dt = meta.dt[file.exists(meta_json),]
    if(!is.null(filter_patients)) {
        meta.dt = meta.dt[patient_id %in% filter_patients,]
    }
    jsons.lst = lapply(1:nrow(meta.dt), function(x) {
        json.dt = jsonlite::read_json(meta.dt$meta_json[x],simplifyVector = TRUE)
    })
    jsons.dt = rbindlist(jsons.lst, fill = TRUE)
                                        #snv distribution json
    snv.dt = jsons.dt[,.(pair, snv_count,tumor_type_final)] %>% setnames(.,c("pair","value","tumor_type_final_mod"))
                                        #sv distribution json
    sv.dt = jsons.dt[,.(pair, sv_count)] %>% setnames(.,c("pair","value"))
    sv.dt[,id := 1:.N]
    sv.dt = sv.dt[,.(id,pair,value)]
                                        #loh
    loh.dt = jsons.dt[,.(pair, tumor_type_final,loh_fraction,loh_seglen,loh_total_genome)] %>% setnames(.,c("pair","tumor_type","value","LOH_seg_len","genome_width"))
                                        #ploidy
    ploidy.dt = jsons.dt[,.(pair, tumor_type_final, ploidy, purity)] %>% setnames(.,c("pair","tumor_type_final","value","purity"))
                                        #purity
    purity.dt = jsons.dt[,.(pair, tumor_type_final, ploidy, purity)] %>% setnames(.,c("pair","tumor_type_final","ploidy","value"))
                                        #coverage variance
    cov_var.dt = jsons.dt[,.(pair, tumor_type_final, dlrs)] %>% setnames(.,c("pair","tumor_type_final_mod","value"))
                                        #tmb
    tmb.dt = jsons.dt[,.(pair, tmb, tumor_type_final)] %>% setnames(.,c("pair","tmb","tumor_type_final"))
                                        #writing jsons
    message(paste0("writing jsons to ",common_folder))
    write_json(snv.dt,paste0(common_folder,"/snvCount.json"),pretty = TRUE)
    write_json(sv.dt,paste0(common_folder,"/svCount.json"),pretty = TRUE)
    write_json(loh.dt,paste0(common_folder,"/lohFraction.json"),pretty = TRUE)
    write_json(ploidy.dt,paste0(common_folder,"/ploidy.json"),pretty = TRUE)
    write_json(purity.dt,paste0(common_folder,"/purity.json"),pretty = TRUE)
    write_json(cov_var.dt,paste0(common_folder,"/coverageVariance.json"),pretty = TRUE)
    write_json(tmb.dt,paste0(common_folder,"/tmb.json"),pretty = TRUE)
}




#' @name bw_temp
#' @title bw_temp
#' @description
#'
#' function to create data.table for bigwig files for pgvdb
#' 
#' @param patient_id patient.id to be added to pgvdb
#' @param order optional entry if you order plots with a column order
#' @param x path to a granges rds or a granges object
#' @param ref reference to use for pgvdb
#' @param chart_type defaultChartType for pgvdb, default is area, can also be scatterplot
#' @param visible TRUE/FALSE whether the plot is hidden or showing in pgv
#' @param title title of the plot in pgvdb
#' @param type bigwig, should not change for this plot type
#' @param field field to plot as the y value, default is foreground
#' @param overwrite TRUE/FALSE to overwrite an existing bw
#' @return NULL
#' @export
#' @author Stanley Clarke

bw_temp = function(patient_id = NA,order = NA, x = list(NA), ref = NA, chart_type = "area", visible = TRUE, title = NA, type = "bigwig", field = "foreground", overwrite = FALSE) {
    dt1 = data.table(patient.id = patient_id,
                     visible = visible,
                     x = x,
                     type = type,
                     field = field,
                     ref = ref,
                     title = title,
                     order = order,
                     defaultChartType = chart_type,
                     overwrite = overwrite
                     )
    return(dt1)
}


#' @name arrow_temp
#' @title arrow_temp
#' @description
#'
#' function to create data.table for arrow files for pgvdb
#' 
#' @param patient_id patient.id to be added to pgvdb
#' @param order optional entry if you order plots with a column order
#' @param x path to a granges rds or a granges object
#' @param ref reference to use for pgvdb
#' @param chart_type defaultChartType for pgvdb, default is scatterplot, can also be area
#' @param visible TRUE/FALSE whether the plot is hidden or showing in pgv
#' @param title title of the plot in pgvdb
#' @param type scatterplot, should not change for this plot type
#' @param field field to plot as the y value, default is foreground
#' @param overwrite TRUE/FALSE to overwrite an existing arrow
#' @return NULL
#' @export
#' @author Stanley Clarke

arrow_temp = function(patient_id = NA, order = NA, x = list(NA), ref = NA, chart_type = "scatterplot", visible = TRUE, title = NA, type = "scatterplot", field = "foreground", overwrite = FALSE) {
    dt1 = data.table(patient.id = patient_id,
                     visible = visible,
                     x = x,
                     type = type,
                     field = field,
                     ref = ref,
                     title = title,
                     order = order,
                     defaultChartType = chart_type,
                     overwrite = overwrite
                     )
    return(dt1)
}

#' @name genome_temp
#' @title genome_temp
#' @description
#'
#' function to create data.table for genome graphs for pgvdb
#' 
#' @param patient_id patient.id to be added to pgvdb
#' @param order optional entry if you order plots with a column order
#' @param x path to a JaBbA ggraph object or object itself as a list
#' @param ref reference to use for pgvdb
#' @param visible TRUE/FALSE whether the plot is hidden or showing in pgv
#' @param title title of the plot in pgvdb
#' @param type genome, can be changed to allelic to render allelic graphs
#' @param annotation default is list of SVs, make null if no annotations present in object
#' @param overwrite TRUE/FALSE to overwrite an existing genome json
#' @return NULL
#' @export
#' @author Stanley Clarke

genome_temp = function(patient_id = NA, order = NA, x = list(NA), ref = NA, type = "genome", visible = TRUE, title = NA, annotation = list(c('bfb','chromoplexy','chromothripsis','del','dm','dup','pyrgo','rigma','simple','tic','tyfonas')), overwrite = FALSE) {
                                        #use type = allelic to make a color a genome graph
    dt1 = data.table(patient.id = patient_id,
                     type = type,
                     visible = visible,
                     title = title,
                     x = x,
                     ref = ref,
                     order = order,
                     annotation = annotation,
                     overwrite = overwrite
                     )
    return(dt1)
}

#' @name walks_temp
#' @title walks_temp
#' @description
#'
#' function to create data.table for walks for pgvdb
#' 
#' @param patient_id patient.id to be added to pgvdb
#' @param order optional entry if you order plots with a column order
#' @param x walks object as a list
#' @param ref reference to use for pgvdb
#' @param visible TRUE/FALSE whether the plot is hidden or showing in pgv
#' @param title title of the plot in pgvdb
#' @param type walk, do not change for this plot type
#' @param annotation default is list of SVs, make null if no annotations present in object
#' @param overwrite TRUE/FALSE to overwrite an existing genome json
#' @param tag optional argument, can be binset to override y spacing in pgv
#' @return NULL
#' @export
#' @author Stanley Clarke

walks_temp = function(patient_id = NA, order = NA, x = list(NA), ref = NA, type = "walk", visible = TRUE, title = NA, tag = NA,overwrite = FALSE) {
    dt1 = data.table(patient.id = patient_id, 
                     visible = visible,
                     x = x,
                     order = order,
                     ref = ref,
                     title = title,
                     overwrite = overwrite
                     )
    return(dt1)
}
