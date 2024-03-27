library(rtracklayer)

gr2bw = function(gr, chrom_lengths, score_col_name, output_filepath) {
    #get seqlengths of specified reference as a named vector
    if (!is.null(gr)) {
        bw.gr = gr %>% sortSeqlevels() %>% gr.chr()
    } else {
        stop("GRanges is null")
    }
                                        #check whether seqlengths are longer than seqlengths
    if(any(bw.gr@seqinfo@seqlengths > chrom_lengths[seqnames %in% names(seqlengths(bw.gr))]$end)) {
        warning(paste0("The seqlengths of your GRanges are longer than the seqlengths of ref. \nAre you sure you are using the correct ref? \nprinting and returning seqlengths of ref..."))
        print(chrom_lengths)
        return(chrom_lengths)
    }
    #    if seqlengths are not longer than ref seqlengths then fix it to have the same seqlengths
    bw.gr@seqinfo@seqlengths = chrom_lengths[seqnames %in% bw.gr@seqinfo@seqnames]$end
    bw.gr$score = bw.gr@elementMetadata[[score_col_name]]
    export.bw(object=bw.gr,con=output_filepath)
    print(paste("Bigwig exported to", output_filepath))
}


#testing
## datafile = "~/Projects/higlass_serv/test_input.rds"
## gr1 = readRDS(datafile)
## gr2bw(datafile=datafile,ref="hg38_chr",score.col="foreground")
## gr2bw(datafile=datafile,ref="hg38",score.col="foreground")
## test = gr2bw(datafile=datafile,ref="hg19",score.col="foreground") # fails purposefully because seqlengths of input are longer than seqlengths of hg19; returns seqlengths of specified reference
## test2 = gr2bw(datafile=datafile,ref="hg19_chr",score.col="foreground") # fails purposefully because seqlengths of input are longer than seqlengths of hg19; returns seqlengths of specified reference
## gr2bw(grange = gr1,ref="hg38_chr",score.col="foreground")
## gr2bw(grange = gr1,score.col="foreground") # fails no ref specified



dt2json = function(dt,patient.id,ref,settings,file_name = paste(getwd(),"test.json",sep="/")) {
    #create vector of seqlengths
    settings_data <- jsonlite::fromJSON(settings)
    chrom_lengths <- as.data.table(settings_data$coordinates$sets[[ref]])[,.(chromosome,startPoint,endPoint)]
    colnames(chrom_lengths) = c("seqnames","start","end")
    chrom_lengths[!grepl("chr",seqnames), seqnames := paste0("chr",seqnames)] # weird fix because hg38_chr does not have chr on Y and M
#convert to ggraph and create json
    gr1 = dt2gr(dt[order(seqnames,start),]) %>% sortSeqlevels() %>% gr.chr()
    if(any(gr1@seqinfo@seqlengths > chrom_lengths[seqnames %in% names(seqlengths(gr1))]$end)) {
        stop(paste("the seqlengths of your granges has ranges that are not contained in the seqlengths of",ref))
    }
    jab = gG(nodes=gr1)
#    jab = gGnome::refresh(jab)
    settings_y = list(y_axis = list(title = "copy number",
                                  visible = TRUE))
    node.json = gr2dt(jab$nodes$gr[, "snode.id"])[, .(chromosome = seqnames, startPoint = start, endPoint = end, iid = snode.id, y = 1,title=snode.id,type="interval",strand="*")]
    gg.js = list(intervals = node.json, connections = data.table())
    gg.js = c(list(settings = settings_y), gg.js)
    jsonlite::write_json(gg.js, file_name,
                         pretty=TRUE, auto_unbox=TRUE, digits=4)
    return(file_name)
}

#hits_medbin = readRDS("~/projects/ATAC_HICHIP/hits_per_norm.depth.med.all.score_offset.dt_8_25_23.rds")[1:1000,]
#dt = dt2json(hits_medbin,settings=settings,ref="hg38")
