find_project_root <- function(start = getwd()) {
  candidate <- normalizePath(start, mustWork = TRUE)

  repeat {
    if (file.exists(file.path(candidate, "data", "ponisio2014dataset.csv"))) {
      return(candidate)
    }

    parent <- dirname(candidate)
    if (identical(parent, candidate)) {
      stop("Could not find the project root containing data/ponisio2014dataset.csv.")
    }
    candidate <- parent
  }
}

prepare_study_level_data <- function(project_root = find_project_root()) {
  datfull <- read.csv(
    file.path(project_root, "data", "ponisio2014dataset.csv"),
    header = TRUE,
    stringsAsFactors = FALSE
  )

  datfull$study_ID <- as.numeric(as.factor(paste(
    datfull$Author,
    datfull$Year,
    datfull$Journal,
    sep = "_"
  )))
  datfull$effect_size_ID <- seq_len(nrow(datfull))

  keep <- datfull$Org.fertilizer.type %in% c("plant", "animal") &
    datfull$Crop.type == "cereals"
  dat <- datfull[keep, , drop = FALSE]
  dat$fertilizer <- factor(
    dat$Org.fertilizer.type,
    levels = c("animal", "plant")
  )

  dat <- metafor::escalc(
    measure = "CVR",
    m1i = Mean.Org,
    sd1i = SD.Org,
    n1i = N.Org,
    m2i = Mean.conv,
    sd2i = SD.conv,
    n2i = N.Conv,
    data = dat,
    var.names = c("lnCVR", "var.lnCVR"),
    correct = TRUE
  )

  group_effects <- table(dat$fertilizer)
  study_fertilizer_count <- tapply(
    dat$fertilizer,
    dat$study_ID,
    function(x) length(unique(x))
  )
  mixed_studies <- as.integer(names(study_fertilizer_count)[
    study_fertilizer_count > 1L
  ])
  mixed_effects <- dat$study_ID %in% mixed_studies
  matched_data <- dat[!mixed_effects, , drop = FALSE]

  stopifnot(
    nrow(dat) == 318L,
    length(unique(dat$study_ID)) == 36L,
    unname(group_effects["animal"]) == 134L,
    unname(group_effects["plant"]) == 184L,
    length(mixed_studies) == 6L,
    sum(mixed_effects) == 86L,
    nrow(matched_data) == 232L,
    length(unique(matched_data$study_ID)) == 30L
  )

  list(
    data = dat,
    matched_data = matched_data,
    mixed_studies = mixed_studies,
    counts = list(
      effects = nrow(dat),
      studies = length(unique(dat$study_ID)),
      animal_effects = unname(group_effects["animal"]),
      plant_effects = unname(group_effects["plant"]),
      mixed_studies = length(mixed_studies),
      mixed_effects = sum(mixed_effects),
      matched_effects = nrow(matched_data),
      matched_studies = length(unique(matched_data$study_ID))
    )
  )
}
