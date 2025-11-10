# app.R
library(shiny)
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(plotly)
library(R.utils)
source("code/analysis_variables.R")

# --- UI ---
ui <- navbarPage(
    title = "MutRisk Gene Explorer",

    # --- Tab 1: Gene Explorer (your original app) ---
    tabPanel("Number of cells with driver mutations",
             sidebarLayout(
                 sidebarPanel(
                     selectInput("tissue", "Tissue", choices = c("colon", "blood", "lung"), selected = "colon"),
                     selectInput("gene_name", "Gene", choices = c("APC", "KRAS")),
                     checkboxInput("exposed_conditions", "Include exposed tissues", value = FALSE),
                     actionButton("run", "Compute", class = "btn-primary")
                 ),
                 mainPanel(
                     plotlyOutput("gene_plots", height = 500, width = 700)
                 )
             )
    ),

    # --- Tab 2 ---
    tabPanel("Driver muts/Gene",
             h3("Position-specific driver counts"),
             p("You can add more inputs, tables, or plots here.")
    ),

    # --- Tab 3 ---
    tabPanel("About",
             h3("About"),
             p("Developed by Axel Rosendahl Huber at the BBGLab Barcelona", tags$a(href = "https://bbglab.irbbarcelona.org/"))
    )
)

# --- Server ---
server <- function(input, output, session) {

    boostDM_list = list(
        "colon" = fread("processed_data/colon_boostDM_cancer.txt.gz"),
        "lung" = fread("processed_data/lung_boostDM_cancer.txt.gz"),
        "blood" = fread("processed_data/CH_boostDM_cancer.txt.gz")
    )

    tissues = c("colon", "lung", "blood")

    # Load metadata
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
        exposed_conditions = input$exposed_conditions

        observe({
            select_tissue <- input$tissue
            choices = tissue_genes[select_tissue]
            updateSelectInput(session, inputId = "gene_name",
                              choices = choices,
                              selected = choices[1])
        })

        select_gene = input$gene_name

        cancer_bDM <- boostDM_list[[select_tissue]]
        metadata = metadata_list[[select_tissue]]
        ratios = ratios_list[[select_tissue]]
        expected_rates = expected_rates_list[[select_tissue]]

        colors = tissue_colors[[select_tissue]]
        ratio_gene = ratios |> filter(gene_name == select_gene)

        ncells_tissue <- tissue_ncells_ci |> filter(tissue %in% select_tissue)
        ncells = ncells_tissue[["mid_estimate"]]

        if (!exposed_conditions) {
            expected_rates = expected_rates |> filter(category %in% c("normal", "non-smoker"))
            metadata = metadata |> filter(category %in% c("normal", "non-smoker"))
        }

        # Gene double vs single
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
}

# --- Launch ---
shinyApp(ui, server)
