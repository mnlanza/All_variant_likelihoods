library(Biostrings)
library(ggplot2)
source("~/Desktop/Relman_lab/All_variant_likelihoods/gyrae_83_all.R")
setwd("~/Desktop/Relman_lab/All_variant_likelihoods")

all_variants <- readDNAStringSet("Double_genes/double_resist_genes.fasta") 
seq_names <- names(all_variants)
"seq_names: double_83_ref double_83_TTT "

# don't reuse has specific file name
get_dfs <- function(){
  dfs <- list()
  for (name in seq_names){
    npy_fn <- paste0("Double_genes/double-resist-genes-v1/output/input_", name, "_logits.npy")
    dfs[[name]] <- initialize_df(npy_fn, name)
  }
  return(dfs)
}

              all_gene_dfs <- get_dfs()
              all_gene_dfs
# returns vector with first - second half of both variants
get_diff_lls_2625 <- function(all_gene_dfs){
  # initialize first half vec
  firsth_ll_vec <- numeric(length(all_gene_dfs))
  firsth_ll_vec
  names(firsth_ll_vec) <- names(all_gene_dfs)
  names(firsth_ll_vec)
  # initialize second half vec
  secondh_ll_vec <- numeric(length(all_gene_dfs))
  secondh_ll_vec
  names(secondh_ll_vec) <- paste0("secondh_", names(all_gene_dfs))
  names(secondh_ll_vec)
  
  for (i in 1:length(all_gene_dfs)) { 
    firsth_ll_vec[i] <- sum(all_gene_dfs[[i]]$log_likelihood[1:2625])
    secondh_ll_vec[i] <- sum(all_gene_dfs[[i]]$log_likelihood[(1+2625):(2625+2625)])
  }
  return(secondh_ll_vec-firsth_ll_vec)
}

dfs <- get_dfs()
get_diff_lls_2625(dfs)


plot_log_likelihood("double_83_ref", dfs$double_83_ref, index_rows  = 200:300, highlighted = 247:249) 

plot_log_likelihood("double_83_ref", dfs$double_83_ref, index_rows  = (200+2625):(300+2625), 
                    highlighted = (2625+247):(2625+249))


plot_log_likelihood("double_83_TTT", dfs$double_83_TTT, index_rows  = 200:300, highlighted = 247:249) 

plot_log_likelihood("double_83_TTT", dfs$double_83_TTT, index_rows  = (200+2625):(300+2625), 
                    highlighted = (2625+247):(2625+249))


plot_log_likelihood("double_83_TTT", dfs$double_83_TTT, index_rows  = (2625):(2625+100), 
                    highlighted = (2600):(2625+100))



##################################
### Making pdfs ###################
##################################
save_as_pdf(plot_log_likelihood("double_83_ref", dfs$double_83_ref, index_rows  = 1:5250, 
            highlighted = 247:249), "2x_83_S1.pdf", out_dir = "Double_genes/My_figures", width = 12, height = 5)

save_as_pdf(plot_log_likelihood("double_83_TTT", dfs$double_83_TTT, index_rows  = 1:5250, 
            highlighted = 247:249) , "2x_83_TTT.pdf", out_dir = "Double_genes/My_figures", width = 12, height = 5)

save_as_pdf(plot_log_likelihood("83_ref", dfs$double_83_ref, index_rows  = 200:300, 
            highlighted = 247:249), "83_ref.pdf", out_dir = "Double_genes/My_figures", width = 12, height = 5)
save_as_pdf(plot_log_likelihood("83_ref_2nd", dfs$double_83_ref, index_rows  = (200+2625):(300+2625), 
            highlighted = (2625+247):(2625+249)), "83_2nd.pdf", out_dir = "Double_genes/My_figures", width = 12, height = 5)




save_as_pdf(plot_log_likelihood("83_TTT_ref", dfs$double_83_TTT, index_rows  = 200:300, 
            highlighted = 247:249), "83_TTT_ref.pdf", out_dir = "Double_genes/My_figures", width = 12, height = 5)
save_as_pdf(plot_log_likelihood("83_TTT_2nd", dfs$double_83_TTT, index_rows  = (2625+200):(2625+300), 
            highlighted = (2625+247):(2625+249)), "83_TTT_2nd.pdf", out_dir = "Double_genes/My_figures", width = 12, height = 5)








save_as_pdf(plot_log_likelihood("trans_83_TTT", dfs$double_83_TTT, index_rows  = (2620):(2600+100),
          highlighted = (2625+247):(2625+249)), "trans_83_TTT.pdf", out_dir = "Double_genes/My_figures", width = 12, height = 5)
save_as_pdf(plot_log_likelihood("trans_83_S1", dfs$double_83_ref, index_rows  = (2620):(2600+100), 
                                highlighted = (2625+247):(2625+249)), "trans_83_S1.pdf", out_dir = "Double_genes/My_figures", width = 12, height = 5)
