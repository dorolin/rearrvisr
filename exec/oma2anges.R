## ------------------------------------------------------------------------
## use OMA output to make ANGES marker input file
##  (this might take a couple of minutes to run)
## ------------------------------------------------------------------------

## input
## -----
## orthologous groups from OMA
ogfile<-"/Users/doro/Dropbox/Drosophila/test/oma/Output/OrthologousGroups.txt"

## only keep subset of species used for ancestral reconstruction
subpfx<-c("ANA","ERE","GRI","MEL","MOJ","PER","PSE","SEC","SIM","VIR","WIL","YAK")
nspec<-length(subpfx)

## genome maps
genomedir<-"/Users/doro/Dropbox/Drosophila/test/anges/input/"
fileext<-rep("_genome.txt",nspec)
genelistcolnames<-c("name","scaff","start","end","strand")
## colnames cannot contain "marker", and require
##  "name","scaff","start","end","strand" (others are ignored)


## ouput
## -----
## markers list for ANGES
outmarkersfile<-"/Users/doro/Dropbox/Drosophila/test/anges/input/Markers"

## modified genome maps
outgenomedir<-"/Users/doro/Dropbox/Drosophila/test/anges/input/"
outfileext<-rep("_markers.txt",nspec)


## processing
## ----------

## avoid having gene positions in output as 3e+05 etc.
op<-options(scipen=999)

## read orthologous groups from OMA
nallspec<-max(count.fields(file=ogfile,sep="\t")) - 1
og.org<-read.table(ogfile,sep="\t",colClasses="character",fill=TRUE,
                   col.names=c("group",paste0("sp",1:nallspec)))

## read genome maps
genelists<-vector("list",nspec)

for(i in 1:nspec){
    genelists[[i]]<-read.table(paste0(genomedir,"/",subpfx[i],fileext[i]),
                               as.is=TRUE)
    colnames(genelists[[i]])<-genelistcolnames
    ## store ortholog IDs below
    genelists[[i]]$marker<-NA
}
names(genelists)<-subpfx

## order orthologs by name
og<-matrix(0,nrow=nrow(og.org),ncol=nspec)
colnames(og)<-subpfx

for(i in 1:nspec){
    mypat<-paste0("^",subpfx[i],":")
    tmp<-matrix("",nrow=nrow(og),ncol=nallspec)
    for(j in 2:ncol(og.org)){
        pos<-grepl(mypat,og.org[,j])
        tmp[pos,j-1]<-og.org[pos,j]
    }
    og[,i]<-apply(tmp,1,function(x) paste(x,collapse=""))
}

## translate the current protein positions (from fasta file)
##  to positions in the genome (check that no NA in match)
om<-matrix(0,nrow=nrow(og),ncol=nspec)
colnames(om)<-subpfx

for(i in 1:nspec){
    genome<-paste0(subpfx[i],":",genelists[[i]]$name)
    a<-match(og[,i],genome)
    a[is.na(a)]<-0
    ## NAs are names in og that are not in genome
    om[,i]<-a
    miss<-which(a==0 & og[,i]!="")
    if(length(miss)!=0){
        warning(paste0("ortholog(s) missing in ",subpfx[i]," genome file (wrong splicing form?):\n  ",paste(og[miss,i],collapse="\n  "),"\ncheck these gene(s) and potentially adjust genome file\n"))
    }
}

## some orthologs might not exist anymore -> exclude
keep<-apply(om,1,function(x) sum(x!=0)>1)
om<-om[keep,,drop=FALSE]
if(nrow(om)<2){
    stop("require at least two orthologs")
}

## make ANGES input file with orthologs
## (runs a little while)
markers<-character()
allmarkers<-character()
prostack<-floor(seq(1,nrow(om),length.out=11)[-1])
proid<-1
cat("\nwriting ANGES input file\n")
for(i in 1:nrow(om)){
    if(i %% prostack[proid] == 0){
        cat(paste0("... ",proid,"0% ...\n"))
        proid<-proid+1
    }
    ## write intermediate steps to improve speed
    if(i %% 100 == 0){
        allmarkers<-c(allmarkers,markers)
        markers<-character()
    }
    markers<-c(markers,paste0(">",i))
    for(j in 1:nspec){
        if(om[i,j]!=0){
            ## get correct genome list
            genome<-genelists[[j]]
            tmp<-genome[om[i,j],]
            tmp2<-paste0(subpfx[j],".",tmp$scaff,":",tmp$start,"-",tmp$end)
            if(tmp$strand == -1){
                tmp2<-paste0(tmp2," -")
            }else if(tmp$strand == "-"){
                tmp2<-paste0(tmp2," -")
            }else if(tmp$strand == 1){
                tmp2<-paste0(tmp2," +")
            }else if(tmp$strand == "+"){
                tmp2<-paste0(tmp2," +")
            }else{
                stop(paste("undefined strand",tmp$strand,"in",tmp$name))
            }
            markers<-c(markers,tmp2)
            ## add ortholog ID to genome
            genelists[[j]]$marker[om[i,j]]<-i
        }
    }
    markers<-c(markers,"")
}
## add last markers
allmarkers<-c(allmarkers,markers)

## re-order columns in genelists and potentially transform strand
for(j in 1:nspec){
    ## re-order
    intercols<-colnames(genelists[[j]])[colnames(genelists[[j]])!="name" & colnames(genelists[[j]])!="marker"]
    genelists[[j]]<-genelists[[j]][c("marker",intercols,"name")]
    ## tansform strand
    genelists[[j]]$strand<-as.character(genelists[[j]]$strand)
    genelists[[j]]$strand<-gsub("^1$","+",genelists[[j]]$strand,perl=TRUE)
    genelists[[j]]$strand<-gsub("^-1$","-",genelists[[j]]$strand,perl=TRUE)
}

## test genomes for unsupported characters in ANGES and
##  overlapping gene positions
for(j in 1:nspec){
    ## get correct genome list
    genome<-genelists[[j]]
    ## only test ANGES genes
    genome<-genome[!is.na(genome$marker),]
    c1<-grep("\\.",genome$scaff)
    c2<-grep(":",genome$scaff)
    c3<-grep("-",genome$scaff)
    c4<-grep("@",genome$scaff)
    c5<-grep("\\|",genome$scaff)
    if(length(c1)!=0){
        warning(paste("unsupported character '.' in scaffold name for",subpfx[j],"\n"))
    }
    if(length(c2)!=0){
        warning(paste("unsupported character ':' in scaffold name for",subpfx[j],"\n"))
    }
    if(length(c3)!=0){
        warning(paste("unsupported character '-' in scaffold name for",subpfx[j],"\n"))
    }
    if(length(c4)!=0){
        warning(paste("unsupported character '@' in scaffold name for",subpfx[j],"\n"))
    }
    if(length(c5)!=0){
        warning(paste("unsupported character '|' in scaffold name for",subpfx[j],"\n"))
    }
    ## test for gene overlaps (ignoring strand)
    genome<-genome[order(genome$scaff,genome$start),]

    overlap<-character()

    i<-1
    pscaff<-genome$scaff[i]
    pend<-genome$end[i]
    pname<-genome$name[i]
    tmpoverlap<-character()

    for(i in 2:nrow(genome)){
        if(genome$scaff[i]==pscaff){
            if(genome$start[i]<=pend){ ## overlap
                if(length(tmpoverlap)==0){
                    tmpoverlap<-paste(pname,genome$name[i],sep="; ")
                }else{
                    tmpoverlap<-paste(tmpoverlap,genome$name[i],sep="; ")
                }
                pend<-max(pend,genome$end[i])
            }else{ ## no overlap
                pend<-genome$end[i]
                pname<-genome$name[i]
                if(length(tmpoverlap)>0){
                    overlap<-c(overlap,tmpoverlap)
                    tmpoverlap<-character()
                }
            }
        }else{ ## next scaffold
            pscaff<-genome$scaff[i]
            pend<-genome$end[i]
            pname<-genome$name[i]
            if(length(tmpoverlap)>0){
                overlap<-c(overlap,tmpoverlap)
                tmpoverlap<-character()
            }
        }
    }
    if(length(tmpoverlap)>0){
        overlap<-c(overlap,tmpoverlap)
    }
    if(length(overlap)>0){
        warning(paste0("overlapping gene positions in ",subpfx[j]," genome file for:\n  ",paste(overlap,collapse="\n  "),"\ncheck these gene positions and adjust accordingly\n"))
    }
}


## write files
## -----------

write.table(allmarkers,file=outmarkersfile,
            quote=FALSE,col.names=FALSE,row.names=FALSE)

for(i in 1:nspec){
    write.table(genelists[[i]],
                file=paste0(outgenomedir,"/",subpfx[i],outfileext[i]),
                quote=FALSE,col.names=TRUE,row.names=FALSE)
}


options(op)

## ------------------------------------------------------------------------

