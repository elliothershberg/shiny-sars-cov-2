library(shiny)
library(tidyverse)
library(r3dmol)
library(JBrowseR)

ui <- fluidPage(
  titlePanel("Sars-CoV-2 Spike Mutations"),
  JBrowseROutput("browserOutput"),
  tags$br(),
  DT::dataTableOutput("mutationsTable"),
  r3dmolOutput("spikeProtein")

)

server <- function(input, output) {
  mutation_data <- read_tsv("spike-mutations.tsv")

  output$mutationsTable <- DT::renderDataTable(
    {
      mutation_data
    },
    selection = "single"
  )

  table_selection <- reactive({
    print(mutation_data[input$mutationsTable_rows_selected, ]$mutation_pos)
    print(mutation_data[input$mutationsTable_rows_selected, ]$mutation_aa)
    mutation_data[input$mutationsTable_rows_selected, ]
  })

  output$spikeProtein <- renderR3dmol(
    spike_model <- r3dmol() %>%
      m_add_model(data = m_fetch_pdb("6vxx"), format = "pdb") %>%
      m_set_style(style = m_style_cartoon()) %>%
      m_add_style(
        style = c(
          m_style_stick(),
          m_style_sphere(scale = 0.3)
        ),
        sel = m_sel(
          resi = table_selection()$mutation_pos
        )
      ) %>%
      m_zoom_to()
  )

  # create the necessary JB2 assembly configuration
  assembly <- assembly(
    "https://jbrowse.org/genomes/sars-cov2/fasta/sars-cov2.fa.gz",
    bgzip = TRUE
  )

  # create configuration for a JB2 GFF FeatureTrack
  annotations_track <- track_feature(
    "https://jbrowse.org/genomes/sars-cov2/sars-cov2-annotations.sorted.gff.gz",
    assembly
  )

  # create the tracks array to pass to browser
  tracks <- tracks(
    annotations_track
  )

  # set up the default session for the browser
  default_session <- default_session(
    assembly,
    c(annotations_track)
  )

  theme <- theme("#5da8a3", "#333")

  location <- reactive({
    if (is.na(table_selection()$chrom_start[1])) {
      "NC_045512.2:1..100"
    } else {
      str_glue(
        "NC_045512.2:",
        "{table_selection()$chrom_start}",
        "..",
        "{table_selection()$chrom_end}"
      )
    }
  })


  # link the UI with the browser widget
  output$browserOutput <- renderJBrowseR(
    JBrowseR(
      "View",
      assembly = assembly,
      tracks = tracks,
      location = location(),
      defaultSession = default_session,
      theme = theme
    )
  )
}

shinyApp(ui = ui, server = server)
