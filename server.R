.libPaths("/local/R/library")

library(ggvis)
library(dplyr)
source("hdf5_utils.R")
if (FALSE) library(RSQLite)

#local
#hdf5_dir <- "/Users/mburger/Documents/mtools/dplot/hdf5/"
#dmz-datasci-srv
hdf5_dir <- "/local/data/hdf5/"

taiga_id <- list(gene.rpkm="013458ff-bb36-4fd7-ac87-dcc04b15174e",
                 demeter="cd8305ab-9ca9-4bac-8d32-d5b0fb499681",
                 demeter.ranks="8bbc1976-0b0f-4e98-a9b0-705afe461989",
                 crispr="539379dd-8588-46d6-9c4a-0bcaf467bbf2",
                 crispr.ranks="1f84eb06-e91c-44de-b05f-7974aabfdfe9",
                 CN="20f80986-b2a9-4ced-948b-3fe0bc0522d7",
                 CTD2="55098e9e-7807-4df6-9c14-bbced64c06bd",
                 in.sets="2256678c-8fe0-4f59-a7d7-a2c7c716bf83",
                 in.sets.cols="cbd33823-2951-4070-adb8-b095dab9a628",
                 dna.missense.inframe="047c1688-1c61-405c-8990-b3bfb9e0e2a3",
                 dna.damage="f9b17f28-0b4c-4e60-89a8-ff9501ab8b2d",
                 dna.hotspot="8a771b60-ef2d-40ff-ad2d-e2a17d42c4b4",
                 rna.missense.inframe="92c22ffa-6204-4808-a767-893011d54f13",
                 rna.damage="f3fd436d-7642-4563-bc9e-28dc32457b7e",
                 rna.hotspot="1d762322-6abb-4cd0-ac23-03d8369492b6",
                 lineage.matrix="ddbb3f5a-62cd-4d0f-ad05-27185cc18772"
                 
  )

mut_combo <- list(dna.mut=paste(taiga_id[["dna.missense.inframe"]],taiga_id[["dna.damage"]],taiga_id[["dna.hotspot"]],sep=","),
                   rna.mut=paste(taiga_id[["rna.missense.inframe"]],taiga_id[["rna.damage"]],taiga_id[["rna.hotspot"]],sep=","))

sample_info <- data.frame(ID=get_values(taiga_id[["in.sets"]],"CCLE_ID",hdf5_dir),lineage=get_values(taiga_id[["in.sets"]],"lineage",hdf5_dir))

axis_vars <- c(
  "Achilles RNAi v2.20 score" = "demeter",
  "Achilles RNAi v2.20 rank" = "demeter.ranks",
  "Achilles CRISPR v3.3.1 score" = "crispr",
  "Achilles CRISPR v3.3.1 rank" = "crispr.ranks",
  "CCLE RNAseq expression" = "gene.rpkm",
  "CCLE CN" = "CN",
  "CTD2 AUC" = "CTD2"
)

default_features <- list(demeter="BRAF",demeter.ranks="BRAF",crispr="BRAF_1_001111",crispr.ranks="BRAF_1_001111",gene.rpkm="BRAF",CN="BRAF",CTD2="nutlin-3")

shinyServer(function(input, output, session) {
   
  x_dataset <- reactive({
    input$x_dataset
  })
    
  observe({
    x_dataset <- input$x_dataset
    x_col <- get_colnames(taiga_id[[x_dataset]],hdf5_dir)
    updateSelectizeInput(session, 'x_gene', choices = x_col, server = TRUE)
  })
  
  observe({
    y_dataset <- input$y_dataset
    y_col <- get_colnames(taiga_id[[y_dataset]],hdf5_dir)
    updateSelectizeInput(session, 'y_gene', choices = y_col, server = TRUE)
  })
  

  
  
  
  
  #output$x_id <- renderUI({
  #  x_cols <- get_colnames(taiga_id[[input$x_dataset]],hdf5_dir)
  #  selectizeInput("x_gene", "Select feature", choices = x_cols, options = list(
  #    placeholder = 'Search',
  #    onInitialize = I('function() { this.setValue(""); }')))
  #})
  
  #output$y_id <- renderUI({
  #  y_cols <- get_colnames(taiga_id[[input$y_dataset]],hdf5_dir)
  #  selectizeInput("y_gene", "Select feature", choices = y_cols, options = list(
  #    placeholder = 'Search',
  #    onInitialize = I('function() { this.setValue(""); }')))
  #})
  
  # Function for generating tooltip text
  x_dataset <- reactive({input$x_dataset})

  # A reactive expression with the ggvis plot
  vis <- reactive({
    
    
    
    m_gene <- "SOX10"
    
    x_dataset <- input$x_dataset
    x_gene <- default_features[[x_dataset]]
    y_dataset <- input$y_dataset
    y_gene <- default_features[[y_dataset]]
    
    if (input$m_dataset == "dna.mut"){
      color_by_datasets <- c("dna.missense.inframe","dna.damage","dna.hotspot")
    } else if (input$m_dataset == "rna.mut"){
      color_by_datasets <- c("rna.missense.inframe","rna.damage","rna.hotspot")
    } else {
      color_by_datasets <- "lineage.matrix"
    }
    
    if ("x_gene" %in% names(input) && nchar(input$x_gene) > 0){
      x_cols <- get_colnames(taiga_id[[x_dataset]],hdf5_dir)
      if (input$x_gene %in% x_cols){
        x_gene <- input$x_gene
      }
    } 
    
    if ("y_gene" %in% names(input) && nchar(input$y_gene) > 0){
      y_cols <- get_colnames(taiga_id[[y_dataset]],hdf5_dir)
      if(input$y_gene %in% y_cols){
        y_gene <- input$y_gene
      }
    }
    
      x <- get_values(taiga_id[[x_dataset]],x_gene,hdf5_dir)
      y <- get_values(taiga_id[[y_dataset]],y_gene,hdf5_dir)
    x_axis_label <- paste0(x_gene," ",get_label(taiga_id[[x_dataset]],hdf5_dir))
    y_axis_label <- paste0(y_gene," ",get_label(taiga_id[[y_dataset]],hdf5_dir))
    
    if (input$filter != "all"){
      
      if (input$filter == "solid" || input$filter == "liquid"){
        indexes <- grep(".+_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE$", names(x))
        liquid_names <- names(x[indexes])
        liquid_samples <- names(x) %in% liquid_names
        if (input$filter == "solid"){
          x <- x[!liquid_samples]
        } else if (input$filter == "liquid"){
          x <- x[liquid_samples]
        }
      } else {
        all_lineages <- get_colnames(taiga_id[["lineage.matrix"]],hdf5_dir)
        if (input$filter %in% all_lineages){
          lin <- get_values(taiga_id[["lineage.matrix"]],input$filter,hdf5_dir)
          lin_samples <- names(lin)[lin == 1]
          x <- x[names(x) %in% lin_samples]
        }
      }
    }
      
      x<-x[! is.na(x)]
      y<-y[! is.na(y)]
      jsamples<-intersect(names(x), names(y))
      xj<-x[jsamples]
      yj<-y[jsamples]
      
      df <- data.frame(ID=jsamples,x_val=xj,y_val=yj,colors=rep("other",length(jsamples)),stringsAsFactors=F)
    
    colors_labs <- list()
    colors_labs[1] <- "other"
    
    if (color_by_datasets == "lineage.matrix" || ("m_gene" %in% names(input) && nchar(input$m_gene) > 0)){
      
      for (i in 1:length(color_by_datasets)){
        m_col <- get_colnames(taiga_id[[color_by_datasets[i]]],hdf5_dir)
        lab <- get_label(taiga_id[[color_by_datasets[i]]],hdf5_dir)
        if (color_by_datasets[1] == "lineage.matrix"){
          m_gene <- input$lin_name
        } else {
          m_gene <- input$m_gene
        }
          
        colors_labs[i+1] <- paste0(m_gene," ",lab)
          
        if (m_gene %in% m_col){
            col <- get_values(taiga_id[[color_by_datasets[i]]],m_gene,hdf5_dir)
            col <- col[names(col) %in% rownames(df)]
            col <- col[!is.na(col)]
            df[names(col[col > 0]),"colors"] <- paste0(m_gene," ",lab)
          }
        
        
      }
    
    
    
    cpalette <- c("#A0A0A0","#A305F2","#F20548","#0B4596","#FA770C","#268533") 
    cpalette <- rev(cpalette[1:length(colors_labs)])
    
      WT_indexes <- grep("other",df$colors)
    if (length(WT_indexes) > 0){
      all_indexes <- 1:nrow(df)
      color_indexes <- all_indexes[-WT_indexes]
      new_order <- c(WT_indexes,color_indexes)
      df <- df[new_order,]
    }
    
      df %>%
        #ggvis(x = xvar, y = yvar) %>%
        ggvis(x = ~x_val, y = ~y_val) %>%
        layer_points(size := 50, size.hover := 200,
          fillOpacity := 1, fillOpacity.hover := 1,
          fill = ~colors, key := ~ID) %>%
          #key := ~ID) %>%
        add_tooltip(point_tooltip, "hover") %>%
        add_axis("x", title = x_axis_label) %>%
        add_axis("y", title = y_axis_label) %>%
        add_legend("fill", title = "Fill", values = rev(unlist(colors_labs))) %>%
        scale_nominal("fill", domain = rev(unlist(colors_labs)),
          range = cpalette) %>%
        set_options(width = 500, height = 500)
    } else {
      df %>%
        #ggvis(x = xvar, y = yvar) %>%
        ggvis(x = ~x_val, y = ~y_val) %>%
        layer_points(size := 50, size.hover := 200,
                     fillOpacity := 1, fillOpacity.hover := 1,
                     fill := "#A0A0A0", key := ~ID) %>%
        #key := ~ID) %>%
        add_tooltip(point_tooltip, "hover") %>%
        add_axis("x", title = x_axis_label) %>%
        add_axis("y", title = y_axis_label) %>%
        #add_legend("fill", title = "Mutations", values = rev(unlist(colors_labs))) %>%
        #scale_nominal("fill", domain = rev(unlist(colors_labs)),
        #              range = cpalette) %>%
        set_options(width = 500, height = 500)
    }
  })

  point_tooltip <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.null(x$ID)) return(NULL)
    
    # Pick out the movie with this ID
    lin <- sample_info[sample_info$ID == x$ID,"lineage"]
    
    
    paste0("<b>", x$ID, "</b><br>"
           ,paste0("lineage: ",lin),"<br>"
           #,paste0("q value: ",cell_line$qvalue),"<br>"
           #,paste0("neg group median: ",cell_line$negative_median),"<br>"
           #,paste0("pos group median: ",cell_line$positive_median))
    )
  }
  
  vis %>% bind_shiny("plot1")
  #vis %>% bind_shiny("plot2")
  output$n_df <- renderUI({
    
    if (input$m_dataset == "dna.mut" || input$m_dataset == "rna.mut"){
      color_by_param <- mut_combo[[input$m_dataset]]
      color_group_param <- input$m_gene
    } else {
      color_by_param <- taiga_id[["lineage.matrix"]]
      color_group_param <- input$m_dataset
    }
    
    #HTML(paste0('<a href="http://52.10.110.160/dplot?x_dataset=',taiga_id[[input$x_dataset]],'&',
    #                                    'y_dataset=',taiga_id[[input$y_dataset]],'&',
    #                                    'lineage_option=',input$filter,'&',
    #                                    'x_gene=',input$x_gene,'&',
    #                                    'y_gene=',input$y_gene,'&',
    #                                    'color_by=',color_by_param,'&',
    #                                    'color_group=',color_group_param,'&',
    #                                    'threshold=15&',
    #                                    'lineages=',taiga_id[["in.sets"]],'&',
    #                                    'lineage_colors=',taiga_id[["in.sets.cols"]],'&',
    #                                    'rule_list=blank">PDF</a>')) 
})
})
