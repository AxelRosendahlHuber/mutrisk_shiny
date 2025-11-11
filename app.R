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

        # --- Load data
        boostDM_list <- list(
            "colon" = fread("processed_data/colon_boostDM_cancer.txt.gz"),
            "lung"   = fread("processed_data/lung_boostDM_cancer.txt.gz"),
            "blood"  = fread("processed_data/CH_boostDM_cancer.txt.gz")
        )
        tissues <- names(boostDM_list)
        tissue_genes <- lapply(boostDM_list, \(x) unique(x$gene_name))

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

        metadata_list  <- lapply(tissues, \(x)
                                 fread(paste0("processed_data/", x, "_metadata.tsv")) |> distinct())
        expected_rates_list <- lapply(tissues, \(x)
                                      fread(paste0("processed_data/", x, "_expected_rates.tsv.gz")))
        ratios_list <-  lapply(tissues, \(x)
                               fread(paste0("processed_data/",  x, "_mut_ratios.tsv.gz")))
        names(ratios_list) = names(expected_rates_list) = names(metadata_list) = tissues

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
            ratio_gene <- ratios |> filter(gene_name == select_gene)

            ncells_tissue <- tissue_ncells_ci |> filter(tissue %in% select_tissue)
            ncells <- ncells_tissue[["mid_estimate"]]

            if (!condition) {
                expected_rates <- expected_rates |> filter(category %in% c("normal", "non-smoker"))
                metadata <- metadata |> filter(category %in% c("normal", "non-smoker"))
            }

            gene_counts_boostdm <- cancer_bDM[gene_name == select_gene & driver == TRUE,
                                              .N, by = c("gene_name", "mut_type", "driver")]

            gene_single_driver <- expected_rates |>
                left_join(ratio_gene, by = "category") |>
                left_join(metadata, by = c("sampleID", "category", "coverage")) |>
                inner_join(gene_counts_boostdm |> filter(driver), by = "mut_type", relationship = "many-to-many") |>
                mutate(across(c(mle, cilow, cihigh), ~ . * N * ratio)) |>
                group_by(category, donor, age, mut_type) |>
                summarize(across(c(mle, cilow, cihigh), mean), .groups = "drop_last") |>
                summarize(across(c(mle, cilow, cihigh), sum))

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

        # --- Load data
        boostDM_list <- list(
            "colon" = fread("processed_data/colon_boostDM_cancer.txt.gz"),
            "lung"   = fread("processed_data/lung_boostDM_cancer.txt.gz"),
            "blood"  = fread("processed_data/CH_boostDM_cancer.txt.gz")
        )
        tissues <- names(boostDM_list)
        tissue_genes <- lapply(boostDM_list, \(x) unique(x$gene_name))

        metadata_list  <- lapply(tissues, \(x) fread(paste0("processed_data/", x, "_metadata.tsv")) |> distinct())
        expected_rates_list <- lapply(tissues, \(x) fread(paste0("processed_data/", x, "_expected_rates.tsv.gz")))
        ratios_list <-  lapply(tissues, \(x) fread(paste0("processed_data/",  x, "_mut_ratios.tsv.gz")))

        names(ratios_list) = names(expected_rates_list) = names(metadata_list) = tissues
        metadata = rbindlist(metadata_list, idcol = "tissue", use.names = TRUE)
        ratios = rbindlist(ratios_list, idcol = "tissue")

        expected_rates = rbindlist(expected_rates_list, idcol = "tissue", use.names = TRUE) |>
            select(-coverage) |>
            left_join(metadata |> select(-coverage, -sensitivity), by = join_by(tissue, sampleID, category))

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
                    ggh4x::facet_grid2(driver ~ ., strip = strip_themed(
                        background_y = elem_list_rect(fill = c("#C03830", "#707071")),
                        text_y = elem_list_text(colour = c("white"), face = "bold")
                    ))
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

