consolidate_figures <- function(source_prefix_map, dest_dir) {

  # Create destination directory if it does not exist
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
    message(paste("Created destination directory:", dest_dir))
  }

  # Iterate through the named list
  for (src_path in names(source_prefix_map)) {
    prefix <- source_prefix_map[[src_path]]

    if (!dir.exists(src_path)) {
      warning(paste("Source directory not found:", src_path))
      next
    }

    # List all files in the source directory
    all_files <- list.files(src_path, full.names = TRUE)

    # Process only files (exclude sub-directories)
    file_list <- all_files[file.info(all_files)$isdir == FALSE]

    for (file_path in file_list) {
      original_name <- basename(file_path)

      # Construct the new filename: [Prefix]_[OriginalName]
      new_name <- paste0(prefix, "_", original_name)
      dest_path <- file.path(dest_dir, new_name)

      # Copy the file to the destination folder
      if (file.copy(file_path, dest_path, overwrite = TRUE)) {
        message(paste("Copied:", new_name))
      } else {
        warning(paste("Failed to copy:", original_name))
      }

    }
  }
}


consolidate_figures(list('simulation/figures' = '04_sim',
                         'diatom/figures' = '05_diatom'),
                    'submission_figures')
