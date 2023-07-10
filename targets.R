library(randomFunctions)
library(devtools)
library(magrittr)
devtools::load_all("/home/cew54/RP/coloc")
source("dirs.rb")
## source("common.R")

## source("~/E/community_functions_v3.R")
## source("~/E/plot_functions.R")
## source("~/E/finemap_functions.R")
## source("~/E/coloc_functions.R")
## source("~/E/common.R")
library(cowplot)
library(parallel)
options(mc.cores=if(interactive()) { 1 } else { 4 })
library(tarchetypes)

args <- getArgs(default=list(block="chr21_block16")) # ETS2
source("ibd-ets2-james-functions.R")

## library(drake)
tar_plan(
    ibd_raw = read_raw("IBD_DeLange_28067908_1-hg38.tsv.gz"),
  psc_raw = read_raw( "PSC_Ji_27992413_1-hg38.tsv.gz"),
  as_raw = read_raw( "ANS_Cortes_23749187_1-hg38.tsv.gz"),
  tak_raw=read_tak(),
  mafld=get_mafld(unique(c(ibd_raw$hm_BP, psc_raw$hm_BP, as_raw$hm_BP))),
  ibd_dataset=make_ibd_dataset(ibd_raw, mafld),
  psc_dataset=make_ibd_dataset(psc_raw, mafld),
  as_dataset=make_ibd_dataset(as_raw, mafld),
  tak_dataset=make_ibd_dataset(tak_raw, mafld),
  s_ibd=runsusie(ibd_dataset),
  ## s_as=runsusie(as_dataset), # susie not happy with AS - multiple ancestries
  s_psc=runsusie(psc_dataset),
  s_tak=runsusie(tak_dataset),
  eqtl_files=get_eqtl_files(),
  eqtl_data=read_eqtl_data(eqtl_files),
  eqtl_lbf=read_eqtl_lbf(eqtl_files),
  eqtl_cs=read_eqtl_cs(eqtl_files),
  eqtl_s=formatsusie_eqtl(eqtl_lbf,eqtl_cs),
  eqtl_dataset=make_eqtl_dataset(eqtl_data),
  final_data=collate_data(c(list(s_ibd, s_psc, NULL, s_tak), eqtl_s),
                          list(ibd_dataset, psc_dataset, as_dataset, tak_dataset, eqtl_dataset)),
  coloc_results=run_coloc(final_data),
  plot=make_plots(final_data,coloc_results)
)

