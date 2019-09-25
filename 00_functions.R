




# General functions -------------------------------------------------------

char <- function(...){
  # Convert non-standard evalutation names (e.g. column names) to a 
  # character vector.
  # 
  # Args:
  #   ... - a sequence of non-standard evaluation names
  #
  # Result:
  #   a character vector
  #
  # Example:
  #   char(a, b, c) --> c("a", "b", "c")
  
  as.vector(map_chr(rlang::quos(...), rlang::quo_name))
}


# Load and save derived tables --------------------------------------------

save_derived <- function(x, filename, compress = "none"){
  # Save a file into the locally defined derived folder (wrapper around
  # `write_rds()` of the readr package)
  #
  # Args:
  #   x - data.table
  #   filename - the name of the rds file (without the extension '.rds')
  #
  # Result:
  #   NULL
  
  path <- file.path(subfolder, der_dir, str_c(filename, ".rds"))
  write_rds(x, path, compress)
}


load_derived <- function(filename, object_name = filename){
  # Read a file from the locally defined derived folder (wrapper around
  # `read_rds()` of the readr package)
  #
  # Args:
  #   filename - the name of the rds file (without the extension '.rds')
  #
  # Result:
  #   the stored object (invisble)
  
  obj <- read_rds(file.path(subfolder, der_dir, str_c(filename, ".rds")))
  assign(object_name, obj, envir = globalenv())
  
  invisible(obj)
}




# Functions to find records of a certain order ----------------------------

first_record <- function(dt){
  # Get the earliest record for each patient in a data.table
  #
  # Args:
  #   dt - data.table with patient id and eventdate
  #
  # Result:
  #   a data.table with one row per patient, which is the earliest event
  
  
  dt[order(patid, eventdate), .SD[1], by = patid]
}


smpl_rnd <- function(dt, n_max = 50000, seed = 1){
  # Take a random sample of patients
  #
  # Parameters
  # ----------
  #  dt : data.table 
  #    source to sample from; must have patid column
  #
  # Return
  # ------
  #  : data.table
  #    the random subsample
  
  if(length(unique(dt$patid)) < n_max){
    return(dt)
  }
  
  set.seed(seed)
  pats <- sample(unique(dt$patid), size = min(n_max, length(unique(dt$patid))))
  dt[data.table(patid = pats), on = "patid"]
}




# Plotting helpers --------------------------------------------------------

# Label comorbidity facets
comorb_lbl <- function(x){
  c("asthma" = "Asthma", 
    "chd" = "CHD", 
    "ckd" = "CKD", 
    "copd" = "COPD",
    "dm" = "Diabetes",
    "hf" = "Heart failure",
    "pad" = "PAD",
    "stroke" = "Stroke",
    "control" = "Non-comorbid controls",
    "none" = "No comorbidity")[x]
}

comorb_labeller <- ggplot2::labeller(comorb = comorb_lbl)



shift_legend <- function(p) {
  # Shift a ggplot legend from the sidelines to an empty facet
  #
  # Heavily adapted (stolen) from https//stackoverflow.com/questions/
  # 54438495/shift-legend-into-empty-facets-of-a-faceted-plot-in-ggplot2
  
  # to grob
  gp <- ggplotGrob(p)
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  
  # establish name of empty panels
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  names <- empty.facet.panels$name
  # example of names:
  #[1] "panel-3-2" "panel-3-3"
  
  # now we just need a simple call to reposition the legend
  lemon::reposition_legend(p, 'center', panel=names)
}

