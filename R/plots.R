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
## create_somatic_json = function(somatic_snv_cn, out_file, pair, pgv_settings, return_table_pgv = FALSE, meta_keep = NULL, y_col = "est_cn_llrm", ref = "hg19") {
##     som.dt = readRDS(somatic_snv_cn)
##     som.dt = som.dt[!is.na(get(y_col)),]
##     som.dt[start == end, end := end +1]
##     som.dt[, strand := NULL]
##     som.dt[variant.p != "",annotation := paste0("Type: ", annotation, "; Gene: ", gene, "; Variant: ",variant.c, "; Protein_variant: ", variant.p, "; VAF: ",vaf)]
##     som.dt[variant.p == "",annotation := paste0("Type: ", annotation, "; Gene: ", gene, "; Variant: ",variant.c, "; VAF: ",vaf)]
##     dt2json_mut(dt = som.dt, ref = ref,settings = pgv_settings, meta_data = meta_keep, y_col = y_col,
##                 file_name = out_file)
##     if(return_table_pgv) {
##         dt.add = data.table(patient.id = pair, type = "genome",visible = TRUE, title = "Copy Number Mutations", source = "mutations.json")
##     }
## }

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
    node.dt = gr2dt(jab$nodes$gr[, c("snode.id","y_value",meta_data)])
    node.json = node.dt[, .(chromosome = seqnames, startPoint = start, endPoint = end, iid = snode.id,title=snode.id,type="interval", y = y_value, annotation = node.dt$annotation)]
    gg.js = list(intervals = node.json, connections = data.table())
    gg.js = c(list(settings = settings_y), gg.js)
    message(paste0("Writing json to ",file_name))
    jsonlite::write_json(gg.js, file_name,
                         pretty=TRUE, auto_unbox=TRUE, digits=4)
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
            res.cn.dt[!is.na(cn) & !is.na(cn.low) & !is.na(cn.high), Variant := paste0("Total CN:",round(cn,digits = 3),"; CN Minor:",round(cn.low,digits = 3),"; CN Major:",round(cn.high,digits = 3))]
            res.cn.dt[!is.na(cn) & is.na(cn.low) & is.na(cn.high), Variant := paste0("Total CN:",round(cn,digits = 3))]
            res.cn.dt[,c("cn", "cn.high", "cn.low", "width", "strand") := NULL] #make null, already added to Variant
            res.final = rbind(res.mut,res.cn.dt)
        } else {
            res.final = res.mut
            res.final[,c("seqnames", "start", "end") := NULL]
        }
        message(paste0("Writing json to ",out_file))
        write_json(res.final, out_file, pretty=TRUE)
        res.final[,sample := pair]
        if(return_table) {
            return(res.final)
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
#' @param make_dir TRUE/FALSE make the directory for the patient sample if it does not exists
#' @return data.table or NULL
#' @export
#' @author Stanley Clarke, Tanubrata Dey, Joel Rosiene

meta_data_json = function(pair, out_file, coverage, jabba_gg, vcf, svaba_somatic_vcf, tumor_type, disease, primary_site, inferred_sex, karyograph, seqnames_loh = c(1:22), seqnames_genome_width = c(1:22,"X","Y"), write_json = TRUE, overwrite = FALSE, return_table = FALSE, make_dir = FALSE) {
    if(!overwrite && write_json == TRUE) {
        if(file.exists(out_file)) {
            print(paste0('Output already exists! - skipping sample ',pair))
            return(NA)
        }
    }
    ## check if directory exists
    ## get folder
    split_file_path = strsplit(out_file, "/")[[1]]
    folder_path = paste0(split_file_path[1:(length(split_file_path)-1)], collapse = "/")
    if(!make_dir) {
        if(!file.exists(folder_path)) {
            print(paste0('Folder does not exist; skipping sample ', pair,". Use make_dir = TRUE to make directory"))
            return(NA)
        }
    }
    if(make_dir) {
        if(!file.exists(folder_path)) {
            cmd = paste0("mkdir -p ", folder_path)
            print(paste0('Making directory ', folder_path))
            system(cmd)
        }
    }
    meta.dt = data.table(pair = pair, tumor_type = tumor_type, tumor_type_final = tumor_type, disease = disease, primary_site = primary_site, inferred_sex = inferred_sex)
    #get derivate log ratio spread
    meta.dt$dlrs = dlrs(readRDS(coverage)$foreground)
    ##Load this case's counts
    ## meta.dt$snv_count = length(read_vcf(vcf))
    vcf.gr = read_vcf(vcf)
    vcf.gr$ALT = NULL # string set slows it down a lot - don't need it here
    meta.dt$snv_count = length(gr.nochr(vcf.gr) %Q% (seqnames %in% seqnames_genome_width))
    ## Count svs, want to count junctions as well as svs
    gg = readRDS(jabba_gg)
    ## cmd = paste0("module unload java && module load java; module load gatk; gatk CountVariants --QUIET true --verbosity ERROR"," -V ",svaba_somatic_vcf)
    ## meta.dt$sv_count = system(paste(cmd, "2>/dev/null"), intern = TRUE)[2] %>% as.integer() #run the command without printing the java command
    ## count just junctions plus loose divided by 2, for sv counts for now
    gg = readRDS(jabba_gg)
    ## meta.dt$junction_count = nrow(gg$junctions$dt[type != "REF",])
    meta.dt$junction_count = nrow(gg$junctions$dt[type != "REF",])
    meta.dt$loose_count = nrow(as.data.table(gg$loose)[terminal == FALSE,])
    meta.dt[,sv_count := (junction_count + (loose_count / 2))]
                                        #get loh
    nodes.dt = gg$nodes$dt
    nodes.dt[, seqnames := gsub("chr","",seqnames)] #strip chr
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
    seqlengths.dt =suppressWarnings(as.data.table(seqinfo(nodes.gr), keep.rownames = "seqnames")) #had to supress, says other arguments ignored
    seqlengths.dt[, seqnames := gsub("chr","",seqnames)] #strip chr
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

create_distributions = function(case_reports_data_folder,common_folder, filter_patients = NULL, write_jsons = TRUE) {
    files.lst = list.files(case_reports_data_folder)
    files.lst = grep("data",files.lst,invert=TRUE, value = TRUE)
    meta.dt = data.table(meta_json = paste0(case_reports_data_folder,files.lst,"/metadata.json"), patient_id = files.lst)
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
    sv.dt = jsons.dt[,.(pair, sv_count, tumor_type_final)] %>% setnames(.,c("pair","value","tumor_type"))
    sv.dt[,id := 1:.N]
    sv.dt = sv.dt[,.(id,pair,value, tumor_type)]
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
    ##temporary fix to make names more consistant
    snv.dt[, tumor_type := tumor_type_final_mod]
    ploidy.dt[, tumor_type := tumor_type_final]
    purity.dt[, tumor_type := tumor_type_final]
    cov_var.dt[, tumor_type := tumor_type_final_mod]
    tmb.dt[, tumor_type := tumor_type_final]
    if(write_jsons == TRUE) {
                                        #writing jsons
        message(paste0("writing jsons to ",common_folder))
        write_json(snv.dt,paste0(common_folder,"/snvCount.json"),pretty = TRUE)
        write_json(sv.dt,paste0(common_folder,"/svCount.json"),pretty = TRUE)
        write_json(loh.dt,paste0(common_folder,"/lohFraction.json"),pretty = TRUE)
        write_json(ploidy.dt,paste0(common_folder,"/ploidy.json"),pretty = TRUE)
        write_json(purity.dt,paste0(common_folder,"/purity.json"),pretty = TRUE)
        write_json(cov_var.dt,paste0(common_folder,"/coverageVariance.json"),pretty = TRUE)
        write_json(tmb.dt,paste0(common_folder,"/tmb.json"),pretty = TRUE)
    } else {
        return(list(snv.dt,sv.dt,loh.dt,ploidy.dt,purity.dt,cov_var.dt,tmb.dt))
    }
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
#' @param max.cn override max cn of 100
#' @param type genome, can be changed to allelic to render allelic graphs
#' @param annotation default is list of SVs, make null if no annotations present in object
#' @param overwrite TRUE/FALSE to overwrite an existing genome json
#' @return NULL
#' @export
#' @author Stanley Clarke

genome_temp = function(patient_id = NA, order = NA, x = list(NA), ref = NA, type = "genome", visible = TRUE, title = NA, max.cn = NULL,annotation = list(c('bfb','chromoplexy','chromothripsis','del','dm','dup','pyrgo','rigma','simple','tic','tyfonas')), overwrite = FALSE) {
                                        #use type = allelic to make a color a genome graph
    dt1 = data.table(patient.id = patient_id,
                     type = type,
                     visible = visible,
                     title = title,
                     x = x,
                     ref = ref,
                     max.cn = max.cn,
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
                     type = type,
                     order = order,
                     ref = ref,
                     title = title,
                     overwrite = overwrite
                     )
    return(dt1)
}

#' @name mutations_temp
#' @title mutations_temp
#' @description
#'
#' function to create data.table for mutations for pgvdb
#' 
#' @param patient_id patient.id to be added to pgvdb
#' @param order optional entry if you order plots with a column order
#' @param x mutations object as a list
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

mutations_temp = function(patient_id = NA, order = NA, x = list(NA), ref = NA, field, type = "mutations", visible = TRUE, title = NA, tag = NA,overwrite = FALSE) {
    dt1 = data.table(patient.id = patient_id, 
                     visible = visible,
                     type = type,
                     x = x,
                     field = field,
                     order = order,
                     ref = ref,
                     title = title,
                     overwrite = overwrite
                     )
    return(dt1)
}

#' @name ppfit_temp
#' @title ppfit_temp
#' @description
#'
#' function to create data.table for ppfit for pgvdb
#' 
#' @param patient_id patient.id to be added to pgvdb
#' @param order optional entry if you order plots with a column order
#' @param x jabba.rds
#' @param ref reference to use for pgvdb
#' @param visible TRUE/FALSE whether the plot is hidden or showing in pgv
#' @param title title of the plot in pgvdb
#' @param type ppfit, do not change for this plot type
#' @param annotation default is NULL, currently won't work for ppfit
#' @param overwrite TRUE/FALSE to overwrite an existing genome json
#' @param tag optional argument, can be binset to override y spacing in pgv
#' @return NULL
#' @export
#' @author Stanley Clarke

ppfit_temp = function(patient_id = NA, order = NA, x = list(NA), ref = NA, type = "ppfit", visible = TRUE, title = NA, annotation = NULL, overwrite = FALSE) {
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

#' @name create_ppfit_json
#' @title create_ppfit_json
#' @description
#'
#' function to create segmentation plots in case reports
#' 
#' @param jabba_rds list object jabba.rds
#' @param out_file location to write json
#' @param write_json TRUE/FALSE whether to write the json
#' @param return_table TRUE/FALSE whether to return the data
#' @param overwrite TRUE/FALSE whether to overwrite existing file
#' @param cores cores for JaBbA:::segstats
#' @return NULL or segstats table
#' @export
#' @author Stanley Clarke, Tanubrata Dey

create_ppfit_json = function(jabba_gg, path_obj, out_file = NULL, write_json = TRUE, overwrite = FALSE, return_table = FALSE, cores = 1) {
    if(!is.null(out_file)) {
        if(!overwrite) {
            if(file.exists(out_file)) {
                print('Output already exists! - skipping sample')
                return(NA)
            }
        }
    }
    #'drcln_cov = readRDS(thisp$decomposed_cov[i])
    ##make this work with complex where the cov file was not an input and with jabba_gg
    ## x = path_obj %>% sniff %>% inputs %>% select(CovFile, maxna) #get coverage that was used for the jabba run
    inputs.dt = path_obj %>% sniff %>% inputs
    if(!any(grepl("CovFile", names(inputs.dt))) & any(grepl("jabba", names(inputs.dt)))) {
        x = inputs.dt$jabba %>% sniff %>% inputs %>% .[,.(CovFile,maxna)]
    } else if (!any(grepl("CovFile", names(inputs.dt))) & any(grepl("jab", names(inputs.dt)))) {
        x = inputs.dt$jab %>% sniff %>% inputs %>% .[,.(CovFile,maxna)]
    } else {
        x = path_obj %>% sniff %>% inputs %>% .[,.(CovFile,maxna)]
    }
    ## coverage = x$CovFile
    cov = readRDS(x$CovFile)
    if ("ratio" %in% names(mcols(cov))) {
        message(paste0("Raw 'cov.rds' was used as input for JaBbA ",path_obj, ", will consider field as 'ratio''\n"))
        field = "ratio"
    } else if ("foreground" %in% names(mcols(cov))) {
        message(paste0("Drycleaned 'drycleaned.cov.rds' was used as input for JaBbA ",path_obj, ", will consider field as 'foreground''\n"))
        field = "foreground"
    }
    if(!(field %in% c("ratio","foreground"))) {
        stop("Cov file is not clear. Ratio nor foreground in the the columns of the coverage file")
    }
    ##need to replace NaN with NA or JaBbA:::segstats breaks
    if(field == "ratio") {
        cov$ratio = gsub("NaN",NA,cov$ratio) %>% as.numeric
    } else if (field == "foreground") {
        cov$foreground = gsub("NaN",NA,cov$foreground) %>% as.numeric
    }
    segstats = JaBbA:::segstats(jabba_gg$nodes$gr,
                                cov,
                                field = field,
                                prior_weight = 1,
                                max.chunk = 1e8,
                                ## subsample = subsample,
                                mc.cores = cores,
                                verbose = FALSE,
                                max.na = as.numeric(x$maxna),
                                lp = FALSE)
    segstats = gr2dt(segstats)
    names(segstats) = gsub("\\.","_",names(segstats))
    ## if (ncol(segstats) == 30) {
    ##     new_names_segstats = c("chromosome","startPoint","endPoint","strand","width","cn",     
    ##                            "start_ix","end_ix","eslack_in","eslack_out","loose","edges_in",
    ##                            "edges_out","tile_id","seg_id","passed","raw_mean","raw_var",   
    ##                            "nbins","nbins_tot","nbins_nafrac","wbins_nafrac","mean","bad", 
    ##                            "max_na","loess_var","tau_sq_post","post_var","var","sd")
    ##     setnames(segstats, old = names(segstats), new = new_names_segstats)
    ## } else if (ncol(segstats) == 29) {
    ##     new_names_segstats = c("chromosome","startPoint","endPoint","strand","width","cn",     
    ##                            "start_ix","end_ix","eslack_in","eslack_out","loose","edges_in",
    ##                            "edges_out","tile_id","seg_id","raw_mean","raw_var",   
    ##                            "nbins","nbins_tot","nbins_nafrac","wbins_nafrac","mean","bad", 
    ##                            "max_na","loess_var","tau_sq_post","post_var","var","sd")
    ##     setnames(segstats, old = names(segstats), new = new_names_segstats)
    ## } else {
    ##     stop(paste0("The expected number of columns are 29 or 30 for segstats. The number of columns for this file are: ", nol(segstats)))
    ## }
    if(write_json) {
        message(paste0("Writing json to ",out_file))
        write_json(segstats, out_file, pretty = TRUE)
    }
    if(return_table) {
        return(segstats)
    }
}


#' @name cov_abs
#' @title cov_abs
#' @description
#'
#' function to run rel2abs on coverage taking either a ggraph or purity & ploidy
#' 
#' @param dryclean_cov dryclean coverage
#' @param jabba_gg optional jabba ggraph. If null, needs purity or ploidy
#' @param purity optional purity value. If null needs ggraph
#' @param ploidy optional ploidy value. If null needs ggraph
#' @param field column in granges to convert with rel2abs
#' @param new_col new column to add to add for the converted rel2abs
#' @return NULL or segstats table
#' @export
#' @author Stanley Clarke
cov2abs = function(dryclean_cov, jabba_gg = NULL, purity = NULL, ploidy = NULL, field = "foreground", new_col = "foregroundabs") {
    cov_gr = readRDS(dryclean_cov)
    if(!is.null(jabba_gg)) {
        gg = readRDS(jabba_gg)
        purity = gg$meta$purity
        ploidy = gg$meta$ploidy
    }
    if(!is.null(purity) && !is.null(ploidy)) {
        purity = purity
        ploidy = ploidy
    }
    mcols(cov_gr)[new_col] = rel2abs(gr = cov_gr,
                                     purity = purity,
                                     ploidy = ploidy,
                                     field = field
                                     )
    return(cov_gr)
}

#' @name cov2arrow_pgv
#' @title cov2arrow_pgv
#' @description
#'
#' function to create arrow with rel2abs from coverage file. Can take a ggraph or specified purity and ploidy to calculate rel2abs
#' 
#' @param patient.id patient name to add to pgv
#' @param dryclean_cov dryclean coverage
#' @param ref reference for pgv
#' @param field column in granges to convert with rel2abs
#' @param jabba_gg optional jabba ggraph. If null, needs purity or ploidy
#' @param purity optional purity value. If null needs ggraph
#' @param ploidy optional ploidy value. If null needs ggraph
#' @param title title for ggraph
#' @param mask optional mask to use after rebinning
#' @param title title for ggraph
#' @param ref reference to use for pgvdb
#' @param title optional title for plot for pgvdb
#' @param seq.fix optional seqlengths vector to fix granges seqlengths
#' @param chart_type default chart type for plot in pgvdb. scatterplot or area
#' @param visible whether the plot is visible in pgv
#' @param field column to do rel2abs on
#' @param new_col new column after rel2abs
#' @param overwrite whether to overwrite the current bigwig
#' @param order optional order if using a column order in your pgvdb object to sort
#' @param binsize size to rebin coverages, default 1e4
#' @return NULL or segstats table
#' @export
#' @author Stanley Clarke
cov2arrow_pgv = function(patient.id, dryclean_cov, jabba_gg = NULL, purity = NULL, ploidy = NULL, mask = NULL, ref, title = NA, seq.fix = NULL, chart_type = "scatterplot", visible = TRUE, field = "foreground", new_col = "foregroundabs", overwrite = FALSE, order = NA, binsize = 1e4) {
    if(!is.null(jabba_gg)) {
        cov_gr = cov2abs(dryclean_cov, jabba_gg, field = field, new_col = new_col)
    } else if(!is.null(purity) && !is.null(ploidy)) {
        cov_gr = cov2abs(dryclean_cov, purity = purity, ploidy = ploidy)
    } else {
        cov_gr = readRDS(dryclean_cov)
    }
    if(!is.null(mask)) {
        cov_gr = gr.val(cov_gr,mask, "mask")
        cov_gr = cov_gr %Q% (is.na(mask))
        cov_gr$mask = NULL
    }
    cov_gr2 = rebin(cov_gr, binsize, field = field)
    if(length(cov_gr2$foregroundabs[cov_gr2$foregroundabs < 0]) > 0) {
        cov_gr2$foregroundabs[cov_gr2$foregroundabs < 0] = 0
    }    
    add.dt = arrow_temp(patient_id = patient.id, ref = ref, field = field, x = list(cov_gr2), title = title, overwrite = overwrite, order = NA, chart_type = chart_type, visible = visible)
    return(add.dt)
}

#' @name cov2bw_pgv
#' @title cov2bw_pgv
#' @description
#'
#' function to create bigwig with rel2abs from coverage file. Can take a ggraph or specified purity and ploidy to calculate rel2abs
#' 
#' @param patient.id patient name to add to pgv
#' @param dryclean_cov dryclean coverage
#' @param ref reference for pgv
#' @param field column in granges to convert with rel2abs
#' @param jabba_gg optional jabba ggraph. If null, needs purity or ploidy
#' @param purity optional purity value. If null needs ggraph
#' @param ploidy optional ploidy value. If null needs ggraph
#' @param mask mask to remove coverages from
#' @param title title for ggraph
#' @param ref reference to use for pgvdb
#' @param title optional title for plot for pgvdb
#' @param seq.fix optional seqlengths vector to fix granges seqlengths
#' @param chart_type default chart type for plot in pgvdb. scatterplot or area
#' @param visible whether the plot is visible in pgv
#' @param field column to do rel2abs on
#' @param new_col new column after rel2abs
#' @param overwrite whether to overwrite the current bigwig
#' @param order optional order if using a column order in your pgvdb object to sort
#' @param mask optional mask to use after rebinning
#' @return NULL or segstats table
#' @export
#' @author Stanley Clarke

## function to not rebin using higlass but mask
cov2bw_pgv = function(patient.id, dryclean_cov, jabba_gg = NULL, purity = NULL, ploidy = NULL, mask = NULL, ref, title = NA, seq.fix = NULL, chart_type = "scatterplot", visible = TRUE, field = "foreground", new_col = "foregroundabs", overwrite = FALSE, order = NA) {
    if(!is.null(jabba_gg)) {
        cov_gr = cov2abs(dryclean_cov, jabba_gg, field = field, new_col = new_col)
    } else if(!is.null(purity) && !is.null(ploidy)) {
        cov_gr = cov2abs(dryclean_cov, purity = purity, ploidy = ploidy)
    } else {
        cov_gr = readRDS(dryclean_cov)
    }
    ## cov_gr = cov2abs(dryclean_cov, jabba_gg, field = field, new_col = new_col)
    if(!is.null(mask)) {
        cov_gr = gr.val(cov_gr,mask, "mask")
        cov_gr = cov_gr %Q% (is.na(mask))
        cov_gr$mask = NULL
    }
    if(length(cov_gr$foregroundabs[cov_gr$foregroundabs < 0]) > 0) {
        cov_gr$foregroundabs[cov_gr$foregroundabs < 0] = 0
    }
    if(!is.null(seq.fix)) {
        ## fix seqlengths to specified seqlengths
        cov_gr = GRanges(as.data.table(cov_gr),seqlengths = seq.fix) %>% trim()
    }
    add.dt = bw_temp(patient_id = patient.id, ref = ref, field = new_col, x = list(cov_gr), title = title, overwrite = overwrite, order = NA, chart_type = chart_type, visible = visible)
    return(add.dt)
}

