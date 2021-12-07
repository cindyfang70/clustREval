library(shiny)
options(shiny.maxRequestSize=300*1024^2) 
ui <- fluidPage(
  fluidRow(
    column(4, 
    fileInput("sce", "Upload SingleCellExperiment Object"),
    textInput("outputPrefix", "Prefix for output file path"),
    "Pipeline Parameters",
    selectInput("filt", label = "Filtering Method", choices = c("filt.default", "filt.lenient", "filt.mad", "filt.pca", "filt.pca2", "filt.stringent", "none"), multiple=TRUE),
    selectInput("norm", label = "Normalization Method", choices =  c("norm.none", "norm.none.scaled", "norm.scnorm.scaled", "norm.scran", "norm.scran.scaled",
                                                                     "norm.sctransform", "norm.scVI", "norm.seurat", "norm.seuratvst"), multiple=TRUE),
    selectInput("sel", label = "Feature Selection Method", choices =  c("sel.expr", "sel.fromField", "sel.vst")),
    numericInput("selnb", "Number of Features to Select", value = 2000, min = 500, max = 10000),
    selectInput("pca", label = "Dimensionality Reduction Method", choices =  c("seurat.pca", "seurat.pca.noweight")),
    selectInput("ndim", "Number of Dimensions for Dimensionality Reduction", choices=c(5, 10, 15, 20, 25), multiple=TRUE),
    selectInput("res", "Clustering Resolution", choices=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2)),
    actionButton("runPipeline", "Run Pipeline", class = "btn-success")
    ),
    column(4, 
           "Gene Set Enrichment Analysis using MSigDB Hallmark Pathways",
           plotOutput("clusts")
           ),
  column(4,
         "Unsupervised Clustering Metrics for Clustering Output"
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
        gseas <- c(gseas, gsea)
      }
      gseas
      })
  output$clusts <- renderPlot({
    gseas <- clusts()
    plotGeneSetEval(gseas[[1]])
  })
}
shinyApp(ui, server)