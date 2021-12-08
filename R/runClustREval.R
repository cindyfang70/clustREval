library(shiny)
options(shiny.maxRequestSize=300*1024^2) 
ui <- fluidPage(
  fluidRow(
    column(4, 
    fileInput("sce", "Upload SingleCellExperiment Object"),
    textInput("outputPrefix", "Prefix for output file path"),
    "Pipeline Parameters",
    selectInput("filt", label = "Filtering Method", choices = c("filt.default", "filt.lenient", "filt.mad", "filt.pca", "filt.pca2", "filt.stringent", "none"), multiple=TRUE),
    selectInput("norm", label = "Normalization Method", choices =  c("norm.scran",
                                                                     "norm.sctransform", "norm.seurat"), multiple=TRUE),
    selectInput("sel", label = "Feature Selection Method", choices =  c("sel.fromField", "sel.vst")),
    numericInput("selnb", "Number of Features to Select", value = 2000, min = 500, max = 10000),
    selectInput("pca", label = "Dimensionality Reduction Method", choices =  c("seurat.pca", "seurat.pca.noweight")),
    selectInput("ndim", "Number of Dimensions for Dimensionality Reduction", choices=c(5, 10, 15, 20, 25), multiple=TRUE),
    selectInput("res", "Clustering Resolution", choices=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2), multiple=TRUE),
    actionButton("runPipeline", "Run Pipeline", class = "btn-success")
    ),
    column(8, 
          "Gene Set Enrichment Analysis using MSigDB Hallmark Pathways",
          uiOutput("pipelineControls"),
          uiOutput("gseaControls"),
          actionButton("runGsea", "Plot Enrichment Results", class="btn-success"),
          plotOutput("clusts"),
         "Unsupervised Clustering Metrics for Clustering Output",
         actionButton("computeMetrics", "Compute Metrics", class="btn-success"),
         tableOutput("metrics")
         )
  )
)


server <- function(input, output, session) {
  clusts <- reactive({
      input$runPipeline
      req(input$sce)
      clusts <- runPipelineCombs(sce=readRDS(input$sce$datapath), outputPrefix=input$outputPrefix, filt=input$filt, norm=input$norm, sel=input$sel, selnb=input$selnb, dr=input$pca, dims=as.double(input$ndim), res=as.double(input$res))
      GMTPath <- system.file("extdata", "h.all.v7.4.symbols.gmt", package = "clustREval")
      gseas <- c()
      for (i in 1:length(clusts)){
        gsea <- geneSetEval(sce=embryo, clusters=clusts[[i]], gmtPathway = GMTPath)
        gseas[[i]] <- gsea
      }
      gseas
      })
  output$pipelineControls <- renderUI({
    gseas <- clusts()
    pipelines <- seq(1, length(gseas))
    print(pipelines)
    selectInput("npipeline", "Choose Pipeline", choices=pipelines)
  })
  output$gseaControls <- renderUI({
    gseas <- clusts()
    req(input$npipeline)
    print(gseas[[as.numeric(input$npipeline)]])
    clusters <- names(gseas[[as.numeric(input$npipeline)]])
    print("clusters:")
    print(clusters)
    selectInput("nclusts", "Choose Cluster", choices=clusters)
  })
  output$clusts <- renderPlot({
    req(input$runGsea)
    req(input$nclusts)
    req(input$npipeline)
    gseas <- clusts()
    print(gseas[[as.numeric(input$npipeline)]])
    print(gseas[[as.numeric(input$npipeline)]][[input$nclusts]])
    plotGeneSetEval(gseas[[as.numeric(input$npipeline)]][[input$nclusts]])
  })
  output$metrics <- renderTable({
    req(input$npipeline)
    req(input$computeMetrics)
    clusts <- runPipelineCombs(sce=readRDS(input$sce$datapath), outputPrefix=input$outputPrefix, filt=input$filt, norm=input$norm, sel=input$sel, selnb=input$selnb, dr=input$pca, dims=as.double(input$ndim), res=as.double(input$res))
    metrics <- computeUnsupervisedMetrics(sce=readRDS(input$sce$datapath), clusters = clusts[[as.numeric(input$npipeline)]])
    print(metrics)
    metrics
  }, digits=10,rownames=TRUE)


}
shinyApp(ui, server)