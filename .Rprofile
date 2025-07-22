options(vsc.use_httpgd = TRUE)
# Use RStudio graphics device only in interactive sessions
if (interactive()) {
    options(device = "RStudioGD")
}

if (file.exists("~/.Renviron")) {
    readRenviron("~/.Renviron")
}

# Set up CRAN and Bioconductor repositories
if (interactive()) {
    options(device = "RStudioGD")
}
#if (!require("BiocManager", quietly = TRUE)) {
#    install.packages("BiocManager")
#}

# Configure language server
options(languageserver.formatting = TRUE)
options(languageserver.diagnostics = TRUE)
options(languageserver.completion = TRUE)
options(languageserver.symbol = TRUE)
options(languageserver.hover = TRUE)

# -- Removed legacy Bioconductor path tweak that called internal
#    BiocManager:::.getBioC_online_sigrepo() (no longer exists)
#    If you need a custom lib path, add it here with valid code.  

# Force graphics to display in Cursor's plot pane
if (interactive() && Sys.getenv("TERM_PROGRAM") == "vscode") {
    options(vsc.plot = TRUE)
    options(vsc.dev = TRUE)
    options(vsc.browser = TRUE)
} 
 