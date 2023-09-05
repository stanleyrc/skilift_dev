#library(rtracklayer)

gr2bw = function(datafile=NULL,grange=NULL,ref,score.col,
                 file_name = paste(getwd(),"test.bw",sep="/")) {
#         file_name = paste0(getwd(),"/",gsub("rds","bw",tstrsplit(datafile,"/")[[length(tstrsplit(datafile,"/"))]]))) { # will save to the same name as your input file but in your current directory with .bw
    print(paste("using",ref,"if this is the wrong reference this command will fail or your coordiantes will be off in PGV"))
#get seqlengths of specified reference as a named vector
    if ((( is.null(datafile) & is.null(grange)) | !is.null(datafile) & !is.null(grange))) { #| datafile != NULL & grange != NULL) {
        stop("to create a bigwig: either an rds of a grange or a grange object must be passed and not both")
    }
    settings_data <- jsonlite::fromJSON(settings)
    chrom_lengths <- as.data.table(settings_data$coordinates$sets[[ref]])[,.(chromosome,startPoint,endPoint)]
    colnames(chrom_lengths) = c("seqnames","start","end")
    chrom_lengths[!grepl("chr",seqnames), seqnames := paste0("chr",seqnames)] # weird fix because hg38_chr does not have chr on Y and M
                                        #creating bigwig
    if (!is.null(datafile)) {
        bw.gr = readRDS(datafile) %>% sortSeqlevels() %>% gr.chr()
    }
    if (!is.null(grange)) {
        bw.gr = grange %>% sortSeqlevels() %>% gr.chr()
    }
#check whether seqlengths are longer than seqlengths
    if(any(bw.gr@seqinfo@seqlengths > chrom_lengths[seqnames %in% names(seqlengths(bw.gr))]$end)) {
#        return(chrom_lengths) ;
        warning(paste0("the seqlengths of your granges is longer than the seqlengths of ",ref,"\nreturning seqlengths of ",ref))
        return(chrom_lengths)
#        return(chrom_lengths) ; stop(paste("the seqlengths of your granges is longer than the seqlengths of ",ref))
        }
#    if seqlengths are not longer than ref seqlengths then fix it to have the same seqlengths
    bw.gr@seqinfo@seqlengths = chrom_lengths[seqnames %in% bw.gr@seqinfo@seqnames]$end
        bw.gr$score = bw.gr@elementMetadata[[score.col]]
        export.bw(object=bw.gr,con=file_name)
        print(paste("bigwig exported to",file_name))
    return(file_name)
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
