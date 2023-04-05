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

args <- getArgs(default=list(block="chr21_block16")) # ETS2

library(drake)
plan <- drake_plan(
  data_raw = read_ibd(),
  mafld=get_mafld(data_raw),
  dataset=make_ibd_dataset(data_raw, mafld),
  s=runsusie(dataset),
  eqtl_files=get_eqtl_files(),
  eqtl_data=read_eqtl_data(eqtl_files),
  eqtl_lbf=read_eqtl_lbf(eqtl_files),
  eqtl_cs=read_eqtl_cs(eqtl_files),
  eqtl_s=formatsusie_eqtl(eqtl_lbf,eqtl_cs),
  eqtl_dataset=make_eqtl_dataset(eqtl_data),
  final_data=collate_data(s, dataset, eqtl_s, eqtl_dataset),
  coloc_results=run_coloc(final_data),
  plot=make_plots(final_data,coloc_results)
)

source("ibd-ets2-james-functions.R")
make(plan)
