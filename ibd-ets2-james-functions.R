
## prepare summary data from IBD
read_raw=function(f) {
  infile=file.path(ROOTDIR, f)
  nm=paste0("zcat ",infile," | head ") %>% fread(cmd=.)
  data=paste0("zcat ",infile," | awk '$3==21 && $4 > 37500000 && $4 < 40000000'") %>%
    fread(cmd=.)
  setnames(data, names(nm))
  if(all(is.na(data$hm_BETA))) # PSC
      data[,hm_BETA:=log(hm_OR)]
  return(data)
}

read_tak=function() {
    infile="/home/cew54/share/Data/GWAS-summary/OrtizFernandez2021-TAK_MetaAnalysis_All5population_Collab.meta.gz"
    if(!file.exists("taktmp/tak.gz-hg38.tsv.gz"))
        system(paste0("mkdir taktmp; zcat ",infile,
                      " | sed 's/A1/REF/' | sed 's/A2/ALT/' |gzip -c > taktmp/tak.gz;",
                      " cd taktmp && ~/share2/GWAS_tools/01-Pipeline/pipeline_v5.3_beta.sh -f tak.gz"))
    nm=paste0("zcat taktmp/tak.gz-hg38.tsv.gz | head  ") %>% fread(cmd=.)
    data=paste0("zcat taktmp/tak.gz-hg38.tsv.gz | awk '$15==21 && $16 > 37500000 && $16 < 40000000'") %>%
        fread(cmd=.)
    setnames(data, names(nm))
    ## table(data$BP38 %in% ibd_raw$hm_BP)
    setnames(data, c("CHR38","BP38","ALT","REF","BETA"),c("hm_CHR","hm_BP","hm_ALT","hm_REF","hm_BETA"))
    ## check alignment
##   chr21_block16 21 37339402 39110976
## chr21_block17 21 39110976 40017600
    DIR="/home/cew54/share2/People/CHRIS/coloc-all"
  refdata16=fread(file.path(DIR,"reference","byblock","chr21_block16.csv"))
  refdata17=fread(file.path(DIR,"reference","byblock","chr21_block17.csv"))
  refdata=rbind(refdata16,refdata17)
  setnames(refdata,c("CHR","BP","REF","ALT","MAF"))
  refdata[,gt:=paste(REF,ALT,sep="/")]
    m=merge(data,refdata,by.x="hm_BP",by.y="BP",suffixes=c(".in",".ref"))
    with(m, table(hm_REF,gt))
   ## hm_REF A/C A/G A/T C/A C/G C/T G/A G/C G/T T/A T/C T/G
   ##   A  63 293  44 151   0   0 674   0   0  96   0   0
   ##   C 123   0   0  63  64 270   0 164   0   0 539   0
   ##   G   0 493   0   0 148   0 230  68  53   0   0 133
   ##   T   0   0  93   0   0 675   0   0 142  45 302  77 
    m[,hm_ALT:=ifelse(hm_REF==REF,ALT,REF)]
    m[,.(hm_CHR,hm_BP,hm_REF,hm_ALT,hm_BETA,SE)]
}

    

## ggplot(data, aes(x=hm_BP, y=-log10(pnorm(-abs(hm_BETA)/SE)*2))) + geom_point() + geom_vline(xintercept = c(38.5e+6,39.5e+6),col="red")


## get LD
get_mafld=function(bp,chr="21") {
    bp  %<>%  unique() 
  ## ss=strsplit(snps,":") %>% do.call("rbind",.) %>%as.data.frame()
  ##   ss$V1 %<>% paste0("chr",.)
  ss=data.table(chr=paste0("chr",chr),
                bp=bp)
  tmp=tempfile()
  tmp.plink=tempfile()
  write.table(ss, file=tmp, quote=FALSE, col.names=FALSE,row.names=FALSE,sep="\t")
  paste("./make_ref_subset.rb ",args$block,tmp,tmp.plink) %>% system()

  ## reset variant identifier in bim file
  file.rename(paste0(tmp.plink,".bim"),paste0(tmp.plink,".bim.orig"))
  paste0("cat ",tmp.plink,".bim.orig | awk 'BEGIN{ OFS=\"\t\" } {print $1,$1 \":\" $4,$3,$4,$5,$6}' > ",tmp.plink,".bim") %>%
    system()

  ## read
  library(annotSnpStats)
  ref16=annot.read.plink(tmp.plink)

  tmp.plink=tempfile()
  newblock="chr21_block17"
  paste("./make_ref_subset.rb ",newblock,tmp,tmp.plink) %>% system()

  ## reset variant identifier in bim file
  file.rename(paste0(tmp.plink,".bim"),paste0(tmp.plink,".bim.orig"))
  paste0("cat ",tmp.plink,".bim.orig | awk 'BEGIN{ OFS=\"\t\" } {print $1,$1 \":\" $4,$3,$4,$5,$6}' > ",tmp.plink,".bim") %>%
    system()

  ## read
  ref17=annot.read.plink(tmp.plink)

  ref=new("SnpMatrix",cbind(ref16@.Data,ref17@.Data))

  message("calculating MAF, LD (can be slow)")
  MAF=col.summary(ref)$MAF; names(MAF)=sub("chr","",colnames(ref))
  drop=which(MAF==0)
  if(length(drop)) {
    MAF=MAF[-drop]
    ref=ref[,-drop]
    snps=names(MAF)
  }
  LD=cor(as(ref,"numeric"))
  ## LD=ld(ref,ref,stat="R",symmetric=TRUE)
  dimnames(LD) %<>% lapply(., function(x) sub("chr","",x))

  ss[,pid:=sub("chr","",paste(chr,bp,sep=":"))]
    message("data pid found in MAFLD data:")
    print(table(ss$pid %in% names(MAF)))
  ss=ss[ pid %in% names(MAF)] # & !is.na(hm_BETA)]
  MAF=MAF[ss$pid]
  LD=LD[ss$pid, ss$pid]
  ref_snps=rbind(ref16@snps,ref17@snps); rownames(ref_snps)=sub("chr","",rownames(ref_snps))
  ref_alleles=with(ref_snps[ss$pid,], paste(allele.1, allele.2, sep="/"))
  return(list(pid=names(MAF),MAF=MAF,LD=LD,ref_alleles=ref_alleles))
}

strReverse <- function(x)
    sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

## check alleles
make_ibd_dataset=function(data, MAFLD) {
  data[,alleles:=paste(hm_ALT,hm_REF,sep="/")]
  data[,pid:=paste(hm_CHR,hm_BP,sep=":")]
  data=data[ pid %in% names(MAFLD$MAF) & !is.na(hm_BETA)]
  ## mfld=MAFLD[ match(data$
  message("allele comparison")
  m= match(data$pid, MAFLD$pid) 
  print(table(data=data$alleles, ref=MAFLD$ref_alleles[m]))
  sw=data$alleles==strReverse(MAFLD$ref_alleles[m])
  if(any(sw))
      data$alleles[sw] = MAFLD$ref_alleles[m][sw]
  all.equal(data$alleles, MAFLD$ref_alleles[m]) # 1 mismatch T/TG vs G/T
  data=data[data$alleles==MAFLD$ref_alleles[ m ], ]
  data=data[!is.na(SE)]
  m= match(data$pid, MAFLD$pid) 
  MAF=MAFLD$MAF[m]
  LD=MAFLD$LD[m,m]

  dataset=with(data,
               list(snp=pid,
                          position=hm_BP,
                          beta=hm_BETA,
                          varbeta=SE^2,
                          MAF=unname(MAF),
                          LD=LD,
                          type="cc",
                          N=12160+13145,
                          s=12160/(12160+13145)))
}

## fm=finemap.signals(dataset,method="cond") # 3 signals, approx quality not good (0.7348, should be 1)
## plot_dataset(dataset, s)

get_eqtl_files=function() {
  raw_files=list.files(file.path(ECATDIR)) %>% grep("monocyte",.,value=TRUE)
  lbf_files=list.files(file.path(ECATDIR,"lbf")) %>% grep("monocyte",.,value=TRUE)
  lbf2raw=sub("snp.txt","all.tsv",lbf_files) %>% gsub("lbf_|_SE","",.)
  lbf2raw %in% raw_files
  list(lbf=lbf_files[ lbf2raw %in% raw_files ],
       raw=lbf2raw[ lbf2raw %in% raw_files ])
}

read_eqtl_data=function(eqtl_files) {
  ## make datasets
  DATA=lapply(file.path(ECATDIR,eqtl_files$raw), function(f) {
    tmp=tempfile()
    message(basename(f))
    nm=paste0("zcat ",f," | awk 'NR==1' ") %>% fread(cmd=.)
    command=paste0("zcat ",f," | awk '$2==21 && $3 > 38500000 && $3 < 39500000' > ",tmp)
    system(command)
    result=fread(tmp)
    unlink(tmp)
    setnames(result, names(nm))
    result[molecular_trait_id %in% c("ENSG00000157557", #"ENSG00000235888",
                                     "ILMN_1720158") ]
  })
  names(DATA)=sub(".all.tsv.gz","",eqtl_files$raw)
  DATA
}

read_eqtl_lbf=function(eqtl_files) {
  LBF=lapply(file.path(ECATDIR, "lbf", eqtl_files$lbf), function(f) {
    tmp=tempfile()
    message(basename(f))
    nm=paste0("zcat ",f," | awk 'NR==1' ") %>% fread(cmd=.)
    command=paste0("zcat ",f," |  awk '$3==21 && $4 > 38500000 && $4 < 39500000' > ",tmp)
    system(command)
    result=fread(tmp)
    unlink(tmp)
    setnames(result, names(nm))
    result[molecular_trait_id=="ENSG00000157557" | molecular_trait_id=="ILMN_1720158" ]
  })
  names(LBF)=sub(".all.tsv.gz","",eqtl_files$raw)
  LBF
}


## check is this two CS? - yes!
## BLUEPRINT_ge_monocyte.all.tsv.gz
## > head(LBF[[1]])
##                 V1                  V2 V3       V4        V5         V6 V7 V8
## 1: ENSG00000157557 chr21_38500041_T_TC 21 38500041 -0.965936 -0.9699809  0  0
## 2: ENSG00000157557  chr21_38500719_T_G 21 38500719 -1.815120 -1.2114195  0  0
## 3: ENSG00000157557  chr21_38500772_G_A 21 38500772 -1.694981 -1.4671507  0  0
## 4: ENSG00000157557  chr21_38501032_T_A 21 38501032 -1.581517 -0.7769799  0  0
## 5: ENSG00000157557  chr21_38501160_G_A 21 38501160 -1.245529 -1.3816373  0  0
## 6: ENSG00000157557  chr21_38502057_C_T 21 38502057 -1.245529 -1.3816373  0  0
##    V9 V10 V11 V12 V13 V14
## 1:  0   0   0   0   0   0
## 2:  0   0   0   0   0   0
## 3:  0   0   0   0   0   0
## 4:  0   0   0   0   0   0
## 5:  0   0   0   0   0   0
## 6:  0   0   0   0   0   0

## readset=function(subindex) {
##   message("reading LBF matrix from ",subindex$lbf_file[1])
##   ranges=paste0("-e '",subindex$st,",",subindex$en,"p'") %>% paste(., collapse=" ")
##   cmd=paste0("zcat ",subindex$lbf_file[1], " | sed -n -e '1p' ",ranges)
##   tmp=fread(cmd=cmd)
##   stmp=split(tmp, tmp$molecular_trait_id)[ subindex$ensg ]
##   names(stmp)=paste0(subindex$shortfile, ".", names(stmp))
##   stmp
## }

## INPUT=split(indices, indices$lbf_file) %>%
##   lapply(., readset) %>%
##   do.call("c", .)
## names(INPUT)=sub(".*.txt.gz.","",names(INPUT))

## remove pairs without snp overlap
  frac_overlap=function(snps1,snps2)
    length(intersect(snps1,snps2))/min(length(snps1),length(snps2))
 
read_eqtl_cs=function(eqtl_files) {
  csfiles=sub(".*lbf_","credible_sets/",eqtl_files$lbf) %>%
    sub(".snp.txt.gz",".purity_filtered.txt.gz",.)
    tmp=tempfile()
  DATA=lapply(file.path(ECATDIR, csfiles), function(f) {
    message(basename(f))
    command=paste0("zcat ",f," | awk 'NR==1' > ",tmp)
    system(command)
    nm=fread(tmp)
    command=paste0("zcat ",f," | awk '$1==\"ENSG00000157557\" || $1==\"ILMN_1720158\"' > ",tmp)
    system(command)
    result=fread(tmp)
    if(is.null(result) || !nrow(result) )
      return(NULL)
    setnames(result, names(nm))
    result[,snp:=paste(chromosome,position,sep=":")]
    cs=with(result, split(snp,cs_index))
    for(i in seq_along(cs))
      names(cs[[i]])=cs[[i]]
    cs
  })
    unlink(tmp)
  names(DATA)=sub(".all.tsv.gz","",eqtl_files$raw)
  DATA
}

formatsusie_eqtl=function(LBF,CS) {
  S=lapply(seq_along(LBF), function(i) {
                                        # find credsets
    lbf=LBF[[i]]
    use=!apply(lbf==0,2,all)
    use=which(use) %>% setdiff(., 1:4)
    if(!length(use))
      return(NULL)
    tmp_susie=list(lbf_variable=t(as.matrix(lbf[,use,with=FALSE])),
                   sets=list(cs=CS[[i]], #list(as.list(structure(seq_along(use), names=paste0("L",1:length(use))))),
                             cs_index=seq_along(use)))
    colnames(tmp_susie$lbf_variable)=paste(lbf$chromosome, # chromosome
                                           lbf$position, # position
                                           sep=":")
      rownames(tmp_susie$lbf_variable)=paste0("L",seq_along(use))
      class(tmp_susie)=c(class(tmp_susie),"susie") # to fool coloc.susie
    tmp_susie
  })
  names(S)=names(LBF)
  S
}

make_eqtl_dataset=function(eqtl_data) {
    ## make snps data for ECAT
    DATA=lapply(eqtl_data, function(d)
        list(snp=paste(d$chromosome,d$position,sep=":"),
             position=d$position,
             beta=d$beta,
             varbeta=d$se^2,
             N=d$an[1],
             type="quant"))
    names(DATA)=names(eqtl_data)
    ## rm dup snps
    for(i in seq_along(DATA)) {
        w=which(duplicated(DATA[[i]]$snp))
        if(length(w)) {
            DATA[[i]]$snp=DATA[[i]]$snp[-w]
            DATA[[i]]$position=DATA[[i]]$position[-w]
            DATA[[i]]$beta=DATA[[i]]$beta[-w]
            DATA[[i]]$varbeta=DATA[[i]]$varbeta[-w]
        }
    }
    DATA
}


## ## check snp overlap
## allsnp=DATA[[1]]$snp
## allecat=lapply(ECATDATA, "[[", "snp") %>% unlist() %>% unique()
## table(allecat %in% allsnp)
## table(allsnp %in% allecat)

run_susie_traitlists=
  function(todo, ...) {
    ret=lapply(1:nrow(todo), function(i) {
        message(i)
      a=run_susie(d1=todo$d1[[i]], d2=todo$d2[[i]], ...)
    }) %>% rbindlist(., fill=TRUE)
  }



run_susie=function(d1,d2,DATA,SUSIE) {
  ## d1=trait idx
  ## d2=trait idx
  ## message("susie coloc: ",d1,"\t",d2)
  dstats1 <- DATA[[d1]] #, list(LD=LD, MAF=MAF))
  dstats2 <- DATA[[d2]] #, list(LD=LD, MAF=MAF))
  susie1=SUSIE[[d1]]
  susie2=SUSIE[[d2]]
  if(frac_overlap(dstats1$snp,dstats2$snp)==0)
    return(NULL)
  ## prep work
  cs1=length(susie1$sets$cs)>0
  cs2=length(susie2$sets$cs)>0
  if(!cs1) {
    if(DATA[[d1]]$type=="quant" & !("sdY" %in% names(DATA[[d1]])))
      DATA[[d1]]$sdY=1
    fm=finemap.abf(DATA[[d1]])
    bf1=structure(fm$lABF., names=fm$snp)
  }
  if(!cs2) {
    if(DATA[[d2]]$type=="quant" & !("sdY" %in% names(DATA[[d2]])))
      DATA[[d2]]$sdY=1
    fm=finemap.abf(DATA[[d2]])
    bf2=structure(fm$lABF., names=fm$snp)
  }
  ## do the coloc
  if(cs1 & cs2) {
    tmp=coloc::coloc.susie(susie1, susie2)$summary
  } else if(!cs1 & !cs2) {
    tmp=coloc::coloc.bf_bf(bf1,bf2)$summary
  } else if(cs1 & !cs2) {
    tmp=coloc::coloc.susie_bf(susie1,bf2)$summary
  } else {
    tmp=coloc::coloc.susie_bf(susie2,bf1)$summary
    ## switch trait order
    allnm=names(copy(tmp))
    setnames(tmp,c("hit1","hit2","idx1","idx2"),c("hit2","hit1","idx2","idx1"))
    tmp=tmp[,allnm,with=FALSE]
  }
  ## annotate
  if(!is.null(tmp)) {
    tmp[,c("trait1","trait2"):=list(d1,d2)]
    if("varbeta" %in% names(dstats1))
      tmp[,hit1.margz:=(dstats1$beta/sqrt(dstats1$varbeta))[match(tmp$hit1,dstats1$snp)]]
    if("varbeta" %in% names(dstats2))
      tmp[,hit2.margz:=(dstats2$beta/sqrt(dstats2$varbeta))[match(tmp$hit2,dstats2$snp)]]
  }
  tmp
}

## run_susie= function(d1,d2,DATA,SUSIE.FM,p1=1e-4,p2=1e-4,p12=5e-6,...) {
##   ## d1=trait name
##   ## d2=trait name
##   ## message("susie coloc: ",d1,"\t",d2)
##   dstats1 <- DATA[[d1]] #, list(LD=LD, MAF=MAF))
##   dstats2 <- DATA[[d2]] #, list(LD=LD, MAF=MAF))
##   if(frac_overlap(dstats1$snp,dstats2$snp)==0)
##     return(NULL)
##   snps_common=intersect(dstats1$snp, dstats2$snp)
##   ## if(ind[[d1]]=="susie")
##     snps_common %<>% intersect(., colnames(SUSIE.FM[[d1]]$lbf_variable))
##   ## if(ind[[d2]]=="susie")
##     snps_common %<>% intersect(., colnames(SUSIE.FM[[d2]]$lbf_variable))
##   ## prep work
##   ## get BF matrix if not susie
##   message(d1,"\t",d2,"\tnsnps = ",length(snps_common))
##   makebf=function(nm) {
##     bf1=SUSIE.FM[[nm]]$lbf_variable[ SUSIE.FM[[nm]]$sets$cs_index, , drop=FALSE]
##     bf1[,snps_common,drop=FALSE]
##   }
##   bf1=makebf(d1)
##   bf2=makebf(d2)
##   if(is.null(bf2))
##     return(NULL)

##   prior_values=list(p1=1e-4, p2=1e-4, p12=5e-6)

##   ## do the coloc
##   tmp=coloc::coloc.bf_bf(bf1, bf2,
##                          p1=prior_values$p1, p2=prior_values$p2, p12=prior_values$p12)$summary
##   ## fix idx of any susie
##   if(ind[[d1]]=="susie")
##     tmp$idx1=SUSIE.FM[[d1]]$sets$cs_index[ tmp$idx1 ]
##   if(ind[[d2]]=="susie")
##     tmp$idx2=SUSIE.FM[[d2]]$sets$cs_index[ tmp$idx2 ]

##   ## annotate
##   if(!is.null(tmp)) {
##     tmp[,c("trait1","trait2"):=list(d1,d2)]
##     tmp[,ind1:=ind[d1]][,ind2:=ind[d2]]
##     if(ind[d1]!="ecat")
##       tmp[,hit1.margz:=(dstats1$beta/sqrt(dstats1$varbeta))[match(tmp$hit1,dstats1$snp)]]
##     if(ind[d2]!="ecat")
##       tmp[,hit2.margz:=(dstats2$beta/sqrt(dstats2$varbeta))[match(tmp$hit2,dstats2$snp)]]
##   }
##   tmp
## }

collate_data=function(susies, datasets,nm) {
    list(ind=c(structure(rep("susie",length(susies)),names=nm[1:length(susies)]),
               structure(rep("ecat",length(eqtl_s)), names=nm[-c(1:length(susies))])),
         DATA=datasets,
         SUSIE.FM=susies)
}

run_coloc=function(susies, datasets) {
  ## concatenate
  message("running coloc...")
    todo <- data.table(d1=1, d2=2:length(susies))
    nm=names(datasets)
    ## expand.grid(d1=c(names(data$ind)[1]),
    ##                   d2=c(names(data$ind)[-1]),
    ##                   stringsAsFactors=FALSE)  %>% as.data.table()
    coloc_results=run_susie_traitlists(todo,DATA=datasets,SUSIE=susies)
    coloc_results[,trait1:=nm[trait1]][,trait2:=nm[trait2]]
  if(nrow(coloc_results)) {
    coloc_results <- coloc_results[!is.na(nsnps)]
    coloc_results <- coloc_results[!is.na(PP.H4.abf)]
  }
  message("...done")
  coloc_results
}

make_plots=function(datasets, susies, coloc_results) {
                                        # everything
    t2=unique(c(coloc_results$trait2, grep("Fairfax|Quach",names(datasets),value=TRUE)))
    png("coloc_all.png",height=9,width=15,units="in",res=300)
    par(mfrow=c(4,5))
    plot_dataset(datasets[[1]], susies[[1]], main="IBD" ,xlim=c(38500000,39500000))
    for(ti in t2)
        plot_dataset(datasets[[ti]], susies[[ti]], main=ti, ,xlim=c(38500000,39500000))
    dev.off()
    print(coloc_results[ PP.H3.abf+PP.H4.abf > 0.5,.(trait2,idx2,PP.H1.abf,PP.H2.abf,PP.H3.abf,PP.H4.abf) ])
    fwrite(coloc_results[ ,.(trait2,idx2,PP.H1.abf,PP.H2.abf,PP.H3.abf,PP.H4.abf) ],
           file="coloc_results.csv")
                                        # zoom
    t2= c("Fairfax_2014_microarray_monocyte_naive", "Fairfax_2014_microarray_monocyte_LPS24")
    png("coloc_zoom.png",height=12,width=8,units="in",res=300)
    par(mfrow=c(5,1))
    for(ti in 1:3) 
        plot_dataset(datasets[[ti]], susies[[ti]], main=names(datasets)[ti],xlim=c(39050000,39150000))
    for(ti in t2)
        plot_dataset(datasets[[ti]], susies[[ti]], main=ti,xlim=c(39050000,39150000))
    dev.off()

}

## coloc_results[order(PP.H4.abf), .(trait1,trait2,PP.H3.abf,PP.H4.abf)]


## plots=vector("list",6)
## for(i in 1:6) {
##   plots[[i]]= if(length(SUSIE.FM[[i]]$sets$cs)) {
##                 plot_dataset(DATA[[i]], SUSIE.FM[[i]])
##               } else {
##                 plot_dataset(DATA[[i]])
##               }
##   a=readLines(n=1)}

## cowplot::plot_grid(plotlist=plots)

show_graph=function() {
    require(igraph)
    a=tar_network(targets_only=TRUE)
    g=graph_from_edgelist(as.matrix(a$edges))
    status=as.data.frame(a$vertices)[ match(names(V(g)), a$vertices$name), "status"]
    V(g)$color= structure(c("grey","lightgreen","orange","red"),names=c("outdated","uptodate","started","errored")) [ status ]
    plot(g)
}
