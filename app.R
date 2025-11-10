# app.R
library(shiny)
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(plotly)
library(R.utils)
source("code/analysis_variables.R")

# --------------------
# Gene Explorer Module
# --------------------
geneExplorerUI <- function(id) {
    ns <- NS(id)
    sidebarLayout(
        sidebarPanel(
            selectInput(ns("tissue"), "Tissue", choices = c("colon", "blood", "lung"), selected = "colon"),
            selectInput(ns("gene_name"), "Gene", choices = c("APC", "KRAS")),
            checkboxInput(ns("exposed_conditions"), "Include exposed tissues", value = FALSE),
            actionButton(ns("run"), "Compute", class = "btn-primary")
        ),
        mainPanel(
            plotlyOutput(ns("gene_plots"), height = 600, width = 800)
        )
    )
}

geneExplorerServer <- function(id) {
    moduleServer(id, function(input, output, session) {

        boostDM_list = list(
            "colon" = fread("processed_data/colon_boostDM_cancer.txt.gz"),
            "lung" = fread("processed_data/lung_boostDM_cancer.txt.gz"),
            "blood" = fread("processed_data/CH_boostDM_cancer.txt.gz")
        )

        tissues = c("colon", "lung", "blood")

        metadata_list  <- lapply(tissues, \(x)
                                 fread(paste0("processed_data/", x, "_metadata.tsv")) |> distinct())
        expected_rates_list <- lapply(tissues, \(x)
                                      fread(paste0("processed_data/", x, "_expected_rates.tsv.gz")))
        ratios_list <-  lapply(tissues, \(x)
                               fread(paste0("processed_data/",  x, "_mut_ratios.tsv.gz")))

        names(ratios_list) = names(expected_rates_list) = names(metadata_list) = tissues
        tissue_genes = lapply(boostDM_list, \(x) x[["gene_name"]] |> unique() )

        plots <- eventReactive(input$run, {

            select_tissue <- input$tissue
            exposed_conditions <- input$exposed_conditions

            observe({
                choices <- tissue_genes[[select_tissue]]
                updateSelectInput(session, "gene_name", choices = choices, selected = choices[1])
            })

            select_gene <- input$gene_name
            cancer_bDM <- boostDM_list[[select_tissue]]
            metadata <- metadata_list[[select_tissue]]
            ratios <- ratios_list[[select_tissue]]
            expected_rates <- expected_rates_list[[select_tissue]]

            colors <- tissue_colors[[select_tissue]]
            ratio_gene <- ratios |> filter(gene_name == select_gene)

            ncells_tissue <- tissue_ncells_ci |> filter(tissue %in% select_tissue)
            ncells <- ncells_tissue[["mid_estimate"]]

            if (!exposed_conditions) {
                expected_rates <- expected_rates |> filter(category %in% c("normal", "non-smoker"))
                metadata <- metadata |> filter(category %in% c("normal", "non-smoker"))
            }

            gene_counts_boostdm <- cancer_bDM[gene_name == select_gene & driver == TRUE, .N, by = c("gene_name", "mut_type",  "driver")]

            gene_single_driver <- expected_rates |>
                left_join(ratio_gene, by = "category") |>
                left_join(metadata, by = c("sampleID", "category", "coverage")) |>
                inner_join(gene_counts_boostdm |> filter(driver), by = "mut_type", relationship = "many-to-many") |>
                mutate(across(c(mle, cilow, cihigh), ~ . * N * ratio)) |>
                group_by(category, donor, age,  mut_type) |>
                summarize(across(c(mle, cilow, cihigh), mean), .groups = "drop_last") |>
                summarize(across(c(mle, cilow, cihigh), sum))

            gene_single_mut_plot <- gene_single_driver |>
                mutate(across(c(mle, cilow, cihigh), ~ . * ncells)) |>
                ggplot(aes(x = age, y = mle, color = category)) +
                geom_pointrange(aes(ymin = cilow, ymax = cihigh)) +
                scale_y_continuous(labels = scales::label_comma()) +
                theme_cowplot() +
                scale_color_manual(values = colors) +
                labs(y = 'number of cells with gene single mutations', x = "Age (years)", title = paste0(select_gene, " Cells with driver mutation"))

            plotly::ggplotly(gene_single_mut_plot) %>% config(displayModeBar = FALSE)
        }, ignoreInit = TRUE)

        output$gene_plots <- renderPlotly(plots())
    })
}

# --------------------
# Tab 1 Module (make barplot of mutations across a range of genes)
# --------------------
sitexplorerUI <- function(id) {
    ns <- NS(id)
    tagList(
        h3("Tab 1 Example"),
        numericInput(ns("num"), "Choose a number:", value = 10, min = 1, max = 100),
        textOutput(ns("result"))
    )
}

sitexplorerServer <- function(id) {
    moduleServer(id, function(input, output, session) {

        # # Load gene_of_interest boostdm
        # boostdm = fread("processed_data/boostdm/boostdm_genie_cosmic/pancancer_boostDM_intersect.txt.gz") |>
        #   mutate(driver = ifelse(driver == TRUE, "driver", "non-driver"))# change names for overview
        #
        # boostDM_list = list(
        #   "colon" = fread("processed_data/colon_boostDM_cancer.txt.gz"),
        #   "lung" = fread("processed_data/lung_boostDM_cancer.txt.gz"),
        #   "blood" = fread("processed_data/CH_boostDM_cancer.txt.gz")
        # )
        #
        # tissues = c("colon", "lung", "blood")
        #
        # metadata_list  <- lapply(tissues, \(x)
        #                          fread(paste0("processed_data/", x, "_metadata.tsv")) |> distinct())
        # expected_rates_list <- lapply(tissues, \(x)
        #                               fread(paste0("processed_data/", x, "_expected_rates.tsv.gz")))
        # ratios_list <-  lapply(tissues, \(x)
        #                        fread(paste0("processed_data/",  x, "_mut_ratios.tsv.gz")))
        #
        #
        #

        output$result <- renderText({
            paste("You chose number:", input$num)
        })
    })
}

# --------------------
# Tab 2 Module (template)
# --------------------
tab2UI <- function(id) {
    ns <- NS(id)
    tagList(
        h3("Tab 2 Example"),
        sliderInput(ns("slider"), "Select a value:", min = 0, max = 50, value = 25),
        textOutput(ns("slider_val"))
    )
}

tab2Server <- function(id) {
    moduleServer(id, function(input, output, session) {
        output$slider_val <- renderText({
            paste("Slider value is:", input$slider)
        })
    })
}

# --------------------
# UI
# --------------------
ui <- navbarPage(
    "MutRisk Explorer",
    tabPanel("Mutated cells/gene during aging", geneExplorerUI("Mutated cells/gene during aging")),
    tabPanel("Mutated cells/position across genes", sitexplorerUI("tab1")),
    tabPanel("About", h3("About"),
             p("Developed by Axel Rosendahl Huber at IRB Barcelona.\n
             For more information visit the", tags$a(href = "https://bbglab.irbbarcelona.org/", target = "_blank", "BBGLab website")))
)

# --------------------
# Server
# --------------------
server <- function(input, output, session) {
    geneExplorerServer("Mutated cells/gene during aging")
    sitexplorerServer("Mutated cells/position across genes")
    tab2Server("tab2") # TODO _ this could be removed if no other tabs are added.
}

# --------------------
# Launch App
# --------------------
shinyApp(ui, server)
