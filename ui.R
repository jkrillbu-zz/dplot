.libPaths("/local/R/library")
library(ggvis)
source("hdf5_utils.R")

#local
#hdf5_dir <- "/Users/mburger/Documents/mtools/dplot/hdf5/"
#dmz-datasci-srv
hdf5_dir <- "/local/data/hdf5/"

axis_vars <- c(
  "Achilles RNAi v2.20 score" = "demeter",
  "Achilles RNAi v2.20 rank" = "demeter.ranks",
  "Achilles CRISPR v3.3.1 score" = "crispr",
  "Achilles CRISPR v3.3.1 rank" = "crispr.ranks",
  "CCLE RNAseq expression" = "gene.rpkm",
  "CCLE CN" = "CN",
  "CTD2 AUC" = "CTD2"
)

lineages <- sort(get_colnames("ddbb3f5a-62cd-4d0f-ad05-27185cc18772",hdf5_dir))
color_vars <- c("DNA mutations" = "dna.mut",
  "RNA mutations" = "rna.mut", "Lineage" = "lineage")
names(lineages) <- lineages


# For dropdown menu
actionLink <- function(inputId, ...) {
  tags$a(href='javascript:void',
         id=inputId,
         class='action-button',
         ...)
}

shinyUI(fluidPage(
  titlePanel("Cancer Dependency Viewer"),
  sidebarLayout(position = "left",
                sidebarPanel(

             h4("X Axis:"),
             #textInput("x_gene", "Gene Symbol",value = "BRAF"),
            selectInput("x_dataset", "Select dataset", axis_vars, selected = "demeter"),
            #conditionalPanel(
              #condition = "input.x_dataset == 'demeter' || input.x_dataset == 'demeter.ranks'",
              selectizeInput(
                'x_gene', 'Select Feature',
                choices = NULL,
                options = list(
                  placeholder = 'Search',
                  onInitialize = I('function() { this.setValue(""); }'))
              ),
            #),
            
            #uiOutput("x_id"),
      
             h4("Y Axis:"),
             #textInput("y_gene", "Gene Symbol",value = "BRAF"),
             selectInput("y_dataset", "Select dataset", axis_vars, selected = "gene.rpkm"),
            selectizeInput(
              'y_gene', 'Select Feature',
              choices = NULL,
              options = list(
                    placeholder = 'Search',
                    onInitialize = I('function() { this.setValue(""); }'))
            ),
            #uiOutput("y_id"),
             h4("Color by:"),
             selectInput("m_dataset", "Annotation Type", color_vars, selected = "rna.mut"),
             conditionalPanel(
               condition = "input.m_dataset != 'lineage' ",
               textInput("m_gene", "Gene Symbol",value = "BRAF")
             ),
             conditionalPanel(
               condition = "input.m_dataset == 'lineage' ",
               selectInput("lin_name", "Select Lineage", lineages, selected = "rhabdoid")
             ),
             h4("Filter:"),
             selectInput("filter", "Cell lines", c("All" = "all",
                                                   "Solid" = "solid",
                                                   "Blood" = "liquid",lineages),
                         selected = "all")
             
           
    ),
    mainPanel(
           ggvisOutput("plot1")
           #,htmlOutput("n_df")
           
           
    )
  )
))