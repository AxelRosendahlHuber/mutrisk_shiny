  # app.R
  library(shiny)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(plotly)
  library(R.utils)
  library(ggh4x)
  source("code/analysis_variables.R")
  source("code/barplot_functions.R")

  # ==================================
  # 1. Global Data Loading
  # ==================================

  # --- Load boostDM data ONCE
  boostDM_list <- list(
      "colon" = fread("processed_data/colon_boostDM_cancer.txt.gz"),
      "lung"  = fread("processed_data/lung_boostDM_cancer.txt.gz"),
      "blood"  = fread("processed_data/CH_boostDM_cancer.txt.gz")
  )
  tissues <- names(boostDM_list)
  tissue_genes <- lapply(boostDM_list, \(x) unique(x$gene_name))

  # --- Load metadata, expected_rates, and ratios ONCE
  metadata_list  <- lapply(tissues, \(x)
                           fread(paste0("processed_data/", x, "_metadata.tsv")) |> distinct())
  expected_rates_list <- lapply(tissues, \(x)
                                fread(paste0("processed_data/", x, "_expected_rates.tsv.gz")))
  ratios_list <- lapply(tissues, \(x)
                        fread(paste0("processed_data/",  x, "_mut_ratios.tsv.gz")))
  names(ratios_list) = names(expected_rates_list) = names(metadata_list) = tissues

  # Combined tables for sitexplorer ONCE
  metadata = data.table::rbindlist(metadata_list, idcol = "tissue", use.names = TRUE)
  ratios = data.table::rbindlist(ratios_list, idcol = "tissue")
  expected_rates = data.table::rbindlist(expected_rates_list, idcol = "tissue", use.names = TRUE) |>
      dplyr::select(-coverage) |>
      dplyr::left_join(metadata |> dplyr::select(-coverage, -sensitivity), by = dplyr::join_by(tissue, sampleID, category))


  # Setkeys for faster data loading
  setkey(ratios, gene_name, category)
  setkey(boostDM_list$colon, gene_name)
  setkey(boostDM_list$lung, gene_name)
  setkey(boostDM_list$blood, gene_name)
  setkey(expected_rates, tissue, category)


  # --------------------
  # Gene Explorer Module
  # --------------------
  geneExplorerUI <- function(id) {
      ns <- NS(id)
      sidebarLayout(
          sidebarPanel(
              selectInput(inputId = ns("tissue"), "Tissue",
                          choices = c("colon", "blood", "lung"), selected = "colon"),
              uiOutput(ns("gene_select_ui")),                # dynamic gene select
              checkboxInput(ns("condition"), "Include exposed conditions", value = FALSE),
              actionButton(ns("run"), "Compute", class = "btn-primary")
          ),
          mainPanel(
              plotlyOutput(ns("gene_plots"), height = 600, width = 800)
          )
      )
  }

  geneExplorerServer <- function(id) {
      moduleServer(id, function(input, output, session) {
          ns <- session$ns

          # --- Dynamic gene select input
          output$gene_select_ui <- renderUI({
              sel_tissue <- ifelse(is.null(input$tissue), "colon", input$tissue)
              choices <- tissue_genes[[sel_tissue]]
              if (is.null(choices) || length(choices) == 0) choices <- "no_genes_available"

              selectInput(
                  inputId = ns("gene_name"),
                  label = "Gene",
                  choices = choices,
                  selected = choices[1],
                  selectize = TRUE
              )
          })

         # --- Plot reactive, updates with button or automatically if needed
          plot_data <- reactive({
              req(input$gene_name)  # wait for dynamic input to exist
              select_tissue <- input$tissue
              condition <- input$condition
              select_gene <- input$gene_name

              cancer_bDM <- boostDM_list[[select_tissue]]
              metadata <- metadata_list[[select_tissue]]
              ratios <- ratios_list[[select_tissue]]
              expected_rates <- expected_rates_list[[select_tissue]]
              colors <- tissue_colors[[select_tissue]]
              ratio_gene <- ratios[gene_name == select_gene,]

              ncells_tissue <- tissue_ncells_ci[tissue_ncells_ci$tissue %in% select_tissue, ]
              ncells <- ncells_tissue[["mid_estimate"]]

              if (!condition) {
                  expected_rates <- expected_rates |> filter(category %in% c("normal", "non-smoker"))
                  metadata <- metadata |> filter(category %in% c("normal", "non-smoker"))
              }

              gene_counts_boostdm <- cancer_bDM[gene_name == select_gene & driver == TRUE,
                                                .N, by = c("gene_name", "mut_type", "driver")]

              # using data table syntax for speed (comments below is dplyr syntax)
              gene_single_driver <- expected_rates[ratio_gene, on = "category"
                ][metadata, on = c("sampleID", "category", "coverage")
                ][gene_counts_boostdm,  on = "mut_type",  nomatch = 0, #0 = Key for Inner Join
                ][j = list(category = category, donor = donor, age = age,
                             mut_type = mut_type, mle = mle * N * ratio,
                           cilow = cilow * N * ratio, cihigh = cihigh * N * ratio)
                 ][, lapply(.SD, mean),.SDcols = c("mle", "cilow", "cihigh"),by = .(category, donor, age, mut_type)][
                      , lapply(.SD, sum),  .SDcols = c("mle", "cilow", "cihigh"),  by = .(category, donor, age)]

              # gene_single_driver <- expected_rates |>
              #     left_join(ratio_gene, by = "category") |>
              #     left_join(metadata, by = c("sampleID", "category", "coverage")) |>
              #     inner_join(gene_counts_boostdm |> filter(driver), by = "mut_type", relationship = "many-to-many") |>
              #     mutate(across(c(mle, cilow, cihigh), ~ . * N * ratio)) |>
              #     group_by(category, donor, age, mut_type) |>
              #     summarize(across(c(mle, cilow, cihigh), mean), .groups = "drop_last") |>
              #     summarize(across(c(mle, cilow, cihigh), sum))

              gene_single_mut_plot <- gene_single_driver |>
                  mutate(across(c(mle, cilow, cihigh), ~ . * ncells)) |>
                  ggplot(aes(x = age, y = mle, color = category)) +
                  geom_pointrange(aes(ymin = cilow, ymax = cihigh)) +
                  scale_y_continuous(labels = scales::label_comma()) +
                  theme_cowplot() +
                  scale_color_manual(values = colors) +
                  labs(
                      y = 'number of cells with gene single mutations',
                      x = "Age (years)",
                      title = paste0(select_gene, " Cells with driver mutation")
                  )

              plotly::ggplotly(gene_single_mut_plot) %>% config(displayModeBar = FALSE)
          })

          # --- Render plot
          output$gene_plots <- renderPlotly({
              # optional: only plot on button click
              input$run
              isolate(plot_data())
          })
      })
  }


  # --------------------
  # Tab 1 Module (make barplot of mutations across a range of genes)
  # --------------------
  sitexplorerUI <- function(id) {
      ns <- NS(id)
      sidebarLayout(
          sidebarPanel(
              selectInput(ns("tissue"), "Tissue", choices = c("colon", "blood", "lung"), selected = "colon"),
              uiOutput(ns("gene_select_ui")),       # dynamic gene select
              uiOutput(ns("condition_select_ui")),  # dynamic condition select
              checkboxInput(ns("split_by_driver"), "Split results by driver", value = FALSE),
              actionButton(ns("run"), "Compute", class = "btn-primary")
          ),
          mainPanel(
              plotlyOutput(ns("gene_plots"), height = 600, width = 900)
          )
      )
  }

  sitexplorerServer <- function(id) {
      moduleServer(id, function(input, output, session) {
          ns <- session$ns

          # --- Dynamic gene select
          output$gene_select_ui <- renderUI({
              sel_tissue <- ifelse(is.null(input$tissue), "colon", input$tissue)
              choices <- tissue_genes[[sel_tissue]]
              if (is.null(choices) || length(choices) == 0) choices <- "no_genes_available"
              selectInput(ns("gene_name"), "Gene", choices = choices, selected = "KRAS", selectize = TRUE)
          })

          # --- Dynamic condition select
          output$condition_select_ui <- renderUI({
              sel_tissue <- ifelse(is.null(input$tissue), "colon", input$tissue)
              choices <- unique(metadata[tissue == sel_tissue][["category"]])
              if (is.null(choices) || length(choices) == 0) choices <- "no_conditions"
              selectInput(ns("condition"), "Condition", choices = choices, selected = choices[1], selectize = TRUE)
          })

          # --- Reactive plot
          plot_data <- reactive({
              req(input$gene_name, input$condition)
              select_tissue <- input$tissue
              select_gene <- input$gene_name
              category_select <- input$condition

              barplot_mutrisk <- make_gene_barplot(
                  expected_rates = expected_rates,
                  boostdm = boostDM_list[[select_tissue]],
                  ratios = ratios,
                  metadata = metadata,
                  tissue_select = select_tissue,
                  gene_of_interest = select_gene,
                  tissue_name = select_tissue,
                  category_select = category_select,
                  lollipop_dots = TRUE
              )

              if (input$split_by_driver) {
                  barplot_mutrisk <- barplot_mutrisk +
                      facet_grid(driver ~ .)
              }

              plotly::ggplotly(barplot_mutrisk) %>% config(displayModeBar = FALSE)
          })

          # --- Render plot
          output$gene_plots <- renderPlotly({
              input$run
              isolate(plot_data())  # plot updates when Compute button clicked
          })
      })
  }


  # --------------------
  # UI
  # --------------------
  ui <- navbarPage(
      title = "MutRisk Explorer",
      tabPanel("Mutated cells per gene during aging", geneExplorerUI("geneExplorer")),
      tabPanel("Mutated cells per position across genes", sitexplorerUI("sitexplorer")),
      tabPanel("About", h3("About"),
               p("Developed by Axel Rosendahl Huber at IRB Barcelona.\n
               For more information visit the", tags$a(href = "https://bbglab.irbbarcelona.org/", target = "_blank", "BBGLab website")))
  )

  # --------------------
  # Server
  # --------------------
  server <- function(input, output, session) {
      geneExplorerServer("geneExplorer")
      sitexplorerServer("sitexplorer")
  }

  # --------------------
  # Launch App
  # --------------------
  shinyApp(ui, server)

