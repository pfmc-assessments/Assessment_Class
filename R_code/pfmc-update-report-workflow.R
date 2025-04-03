# Script to generate stock assessment update for west coast groundfish

# Make sure tinytex is installed
install.packages("tinytex")

# Download asar
remotes::install_github("nmfs-ost/asar")
# Optional: download stockplotr
remotes::install_github("nmfs-ost/stockplotr")

# If that doesn't work, try 
# pak::pak("nmfs-ost/asar")
# pak::pak("nmfs-ost/stockplotr")

# If that doesn't work, try 
# install.packages("asar", repos = c("https://nmfs-ost.r-universe.dev", "https://cloud.r-project.org"))
# install.packages("stockplotr", repos = c("https://nmfs-ost.r-universe.dev", "https://cloud.r-project.org"))

# If already have model results
  # asar::convert_output(
  #   output_file = "Report.sso",
  #   outdir = "~/path/to/report/folder",
  #   model = "SS3",
  #   file_save = TRUE,
  #   savedir = "~/path/to/report/folder",
  #   save_name = "yelloweye_rockfish"
  # )

# There is a function asar::create_template()
# I do not recommend using it.
# Instead, copy the template files from the class repository and edit them
# The asar package created an initial template that I edited.
# I am unsure whether you will need to install asar. (I have it installed)

