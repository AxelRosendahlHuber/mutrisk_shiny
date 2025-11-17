
merge_mutrisk_drivers = function(expected_rates, boostdm, ratios, metadata, gene_of_interest, tissue_select = "colon", tissue_name,
                                 category_select = "normal", cell_probabilities = FALSE,
                                 individual = FALSE, older_individuals = TRUE) {

  older_individuals = metadata |> filter(tissue == tissue_select,
                                         category %in% category_select,
                                         age > 30)
  ratio_gene_tissue = ratios |> filter(gene_name == gene_of_interest &
                                         category %in% category_select,
                                       tissue == tissue_select) |> pull(ratio)
  expected_rates_select = expected_rates[category %in% category_select &
                                           tissue == tissue_select, ] |>
    left_join(older_individuals, by = c("tissue", "sampleID", "category", "donor")) |>
    filter(donor %in% older_individuals$donor)

  ncells_select = tissue_ncells_ci[tissue_ncells_ci$tissue == tissue_select, "mid_estimate"]

  if (cell_probabilities == TRUE) {
    ncells_select = 1
  }

  # group by donor individually
  mutated_rates = expected_rates_select |>
    group_by(donor, mut_type, tissue) |>
    summarize(across(c(mle, cilow, cihigh), mean), .groups = "drop") |>
    mutate(across(c(mle, cilow, cihigh), ~ . * ratio_gene_tissue * ncells_select))

  # modify for specific individual, or for the entiriety
  if (individual == FALSE) {
    print("taking mean of the mutation rates")

    mutated_rates_select = mutated_rates |>
      group_by(mut_type, tissue) |>
      summarize(mle = mean(mle))

    individuals = older_individuals |>
      select(donor, age) |> distinct()

    label = paste(tissue_name, "- average age:", format(mean(individuals$age), digits = 3))


  } else if (individual %in% unique(mutated_rates$donor)) {
    mutated_rates_select = mutated_rates |>
      filter(donor == individual) |>
      select(mut_type, tissue, mle)

    label = paste(category_select, "donor", individual, " age:", older_individuals[donor == individual] |> pull(age))

  }  else if (individual == "all") {
    mutated_rates_select = mutated_rates
    label = "no_label"
  }   else {print("parameter 'individual' must either be FALSE, or donor-id")}

  boostdM_goi = boostdm[gene_name == gene_of_interest, c("mut_type", "position", "driver")]
  expected_gene_muts = boostdM_goi |>
    full_join(mutated_rates_select, relationship = "many-to-many", by = "mut_type")

  return(list(expected_gene_muts = expected_gene_muts, label = label))
}



make_gene_barplot = function(expected_rates, boostdm, ratios, metadata, gene_of_interest,
                             tissue_select = "colon", tissue_name = NULL,
                             category_select = "normal",
                             cell_probabilities = FALSE, individual = FALSE, older_individuals = TRUE,
                             lollipop_dots = FALSE) {

  if (is.null(tissue_name)) {tissue_name = tissue_select}

  mr_drivers = merge_mutrisk_drivers(expected_rates, boostdm, ratios, metadata,
                                     gene_of_interest, tissue_select, tissue_name, category_select, cell_probabilities,
                                     individual, older_individuals)

  expected_gene_muts = mr_drivers$expected_gene_muts
  label = mr_drivers$label

  y_label = "Number of cells with mutation"
  if (cell_probabilities == TRUE) {
    ncells_select = 1
    y_label = "Probability of mutation\n per cell(x10⁻⁶)"
  }

  if (max(expected_gene_muts$position, na.rm = TRUE) > Inf) { # for now set the level to Inf to allow for large genes
    expected_gene_muts = expected_gene_muts |>
      mutate(position = (position - 1) %/% 5 + 1,
             position = position * 5) |>
      group_by(position, tissue, mut_type, driver)  |>
      summarise(mle = sum(mle, na.rm = TRUE), .groups = "drop")
    x_label = "AA position (5AA bins)"
  } else { x_label = "AA position"}

  expected_gene_muts_label = left_join(expected_gene_muts, triplet_match_substmodel)
  expected_gene_muts_label = expected_gene_muts_label |>
    filter(!is.na(driver))

  expected_gene_muts_label = expected_gene_muts_label |>
    mutate(driver = ifelse(driver, "driver mutations", "non-driver\n mutation"))


  # way to make the plot extend both upper and lower axes
  pl = ggplot(expected_gene_muts_label,
              aes(x = position, y = mle)) +
    geom_col(aes(fill = type)) +
    scale_fill_manual(values = COLORS6) +
    theme_cowplot() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = label_comma()) +
    scale_x_continuous(expand = expansion(mult = c(0.01,0.01))) +
    labs(x = x_label, y = y_label, title = gene_of_interest, subtitle = label, fill = NULL) +
    theme(legend.position = "bottom")

  if (cell_probabilities == TRUE) {
    pl = pl + scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = function(x) x * 1e6)
  }

  pl
}
