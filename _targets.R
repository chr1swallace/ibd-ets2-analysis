# _targets.R file
library(targets)
library(tarchetypes)
source("ibd-ets2-james-functions.R")
devtools::load_all("/home/cew54/RP/coloc")
source("dirs.rb")
tar_option_set(packages = c("randomFunctions", "devtools", "magrittr","parallel","cowplot"))
## list(
##   tar_target(file, "data.csv", format = "file"),
##   tar_target(data, get_data(file)),
##   tar_target(model, fit_model(data)),
##   tar_target(plot, plot_model(model, data))
## )

options(mc.cores=if(interactive()) { 1 } else { 4 })
library(randomFunctions)
library(devtools)
library(magrittr)
library(parallel)
library(cowplot)
library(data.table)

args <- list(block="chr21_block16") # ETS2

tar_plan(
    ibd_raw       =read_raw("IBD_DeLange_28067908_1-hg38.tsv.gz"),
    psc_raw       =read_raw( "PSC_Ji_27992413_1-hg38.tsv.gz"),
    as_raw        =read_raw( "ANS_Cortes_23749187_1-hg38.tsv.gz"),
    tak_raw=read_tak(),
    mafld         =get_mafld(unique(c(ibd_raw$hm_BP, psc_raw$hm_BP, as_raw$hm_BP))),
    ibd_dataset   =make_ibd_dataset(ibd_raw, mafld),
    psc_dataset   =make_ibd_dataset(psc_raw, mafld),
    as_dataset    =make_ibd_dataset(as_raw, mafld),
    tak_dataset=make_ibd_dataset(tak_raw, mafld),
    s_ibd         =runsusie(ibd_dataset),
    s_psc         =runsusie(psc_dataset),
    s_tak=runsusie(tak_dataset),
    ## NB s_as doesn't work, b/c multiancestry
    eqtl_files    =get_eqtl_files(),
    eqtl_filt     =filter_eqtl_data(eqtl_files),
    eqtl_data     =read_eqtl_data(eqtl_files,eqtl_filt),
    eqtl_lbf      =read_eqtl_lbf(eqtl_files),
    eqtl_cs       =read_eqtl_cs(eqtl_files), 
    eqtl_s        =formatsusie_eqtl(eqtl_lbf,eqtl_cs),
    eqtl_dataset  =make_eqtl_dataset(eqtl_data),
    susies        =c(list(ibd=s_ibd,psc=s_psc,as=NULL,tak=s_tak),eqtl_s),
    datasets      =c(list(ibd=ibd_dataset,psc=psc_dataset,as=as_dataset,tak=tak_dataset),eqtl_dataset),
    coloc_results =run_coloc(susies, datasets),
    plot          =make_plots(datasets,susies,coloc_results)
)

if(FALSE) {

    tar_load(coloc_results)
    options(width=100)
    coloc_results[,.(nsnps,signif(PP.H3.abf,2),signif(PP.H4.abf,2),idx2,trait2,hit2.margz)]

    tmp=fread("coloc_results.csv")
}
