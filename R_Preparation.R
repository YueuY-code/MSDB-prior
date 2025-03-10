######################################################################
#------------------First: Install the dependent packages-------------#
######################################################################
#Package installation function
install_packages_safe <- function(packages) {
  installed <- c()
  failed <- c()
  options(install.packages.compile.from.source = "no")
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing package:", pkg))
      tryCatch({
        install.packages(pkg)
        installed <- c(installed, pkg)
      }, error = function(e) {
        message(paste("Failed to install package:", pkg, "Error:", e$message))
        failed <- c(failed, pkg)
      })
    } else {
      message(paste("Package", pkg, "is already installed. Skipping..."))
    }
  }
  return(failed )
}
#Package dependent
packages_to_install <- c("SAMprior","base","tidyr","MASS","VGAM","RBesT","dplyr","keras",
                         "overlapping","nnet" ,"knitr","survival","distr","tidyr",
                         "zeallot","EBrmap","doParallel","foreach","BH","cmdstanr")
#installation
install_packages_safe(packages_to_install)
foreach::foreach(i=1:length(packages_to_install))%do%{
  library(packages_to_install[i],character.only = T)}
######################################################################
#------------------Secend: Install cmdstan environment---------------#
######################################################################
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE)
#Unzip the cmdstan package to the specified path,eg:D:cmdstan-2.35.0
#Then set the location of cmdstan in R
set_cmdstan_path("D:cmdstan-2.35.0")
#Compile cmdstan in R
rebuild_cmdstan(cores = 10)