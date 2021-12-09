library(shiny)
library(shinycssloaders)
options(shiny.maxRequestSize=300*1024^2) 
ui <- fluidPage(
    mainPanel(
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
            tabPanel("Specify Pipeline Parameters",
                    fluidRow(
                      column(4, 
                            "Please specify your pipeline parameters and run the pipelines on the next tab.",
                            br(), br(),
                            fileInput("sce", "Upload SingleCellExperiment Object as RDS file"),
                            textInput("outputPrefix", "Prefix for output file path")
                      ),
                      column(4,
                            "Pipeline Parameters",
                            selectInput("filt", label = "Filtering Method", choices = c("filt.default", "filt.lenient", "filt.stringent", "none"), multiple=TRUE),
                            selectInput("norm", label = "Normalization Method", choices =  c("norm.scran",
                                                                                                   "norm.sctransform", "norm.seurat"), multiple=TRUE),
                            selectInput("sel", label = "Feature Selection Method", choices =  c("sel.vst")),
                            numericInput("selnb", "Number of Features to Select", value = 2000, min = 500, max = 10000),
                            selectInput("pca", label = "Dimensionality Reduction Method", choices =  c("seurat.pca", "seurat.pca.noweight")),
                            selectInput("ndim", "Number of Dimensions for Dimensionality Reduction", choices=c(5, 10, 15, 20, 25), multiple=TRUE),
                            selectInput("res", "Clustering Resolution", choices=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2), multiple=TRUE),
                            #actionButton("runPipeline", "Run Pipelines", class = "btn-success")
                                    ),
                           column(4,
                                  "Experiment Summary",
                                  tableOutput("sceSummary")
                                  )
                         )
                         ),
                tabPanel("Run Pipeline and Plot GSEA",
                         fluidRow(
                           column(12,
                         "Gene Set Enrichment Analysis using MSigDB Hallmark Pathways",
                         br(),
                         br(),
                         "Pipelines Run",
                         br(),
                         br(),
                         tableOutput("pipelinesRun"),
                         uiOutput("pipelineControls"),
                         withSpinner(uiOutput(outputId = "gseaControls")),
                         actionButton("runGsea", "Run Pipelines and Plot Enrichment Results", class="btn-success"),
                         plotOutput("clusts") %>% withSpinner(color="#0dc5c1")
                        )
                
                         )
                  ),
                tabPanel("Unsupervised Clustering Metrics", 
                         fluidRow(
                           column(8,
                                  "Unsupervised Clustering Metrics for Clustering Output",
                                  br(), br(),
                                  uiOutput("metricsControls"),
                                  actionButton("computeMetrics", "Compute Metrics", class="btn-success"),
                                  tableOutput("metrics"))%>% withSpinner(color="#0dc5c1")
                                  )
                         )

      )
))


server <- function(input, output, session) {
  gseaRes <- reactive({
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
  clustRes <- reactive({
    req(input$sce)
    clusts <- runPipelineCombs(sce=readRDS(input$sce$datapath), outputPrefix=input$outputPrefix, filt=input$filt, norm=input$norm, sel=input$sel, selnb=input$selnb, dr=input$pca, dims=as.double(input$ndim), res=as.double(input$res))
    clusts
  })
  output$pipelinesRun <- renderTable({
    req(input$sce)
    clusts <- clustRes()
    pips <- pipeComp::parsePipNames(names(clusts))
    rownames(pips) <- seq(1:nrow(pips))
    pips
  }, rownames=TRUE)
  output$pipelineControls <- renderUI({
    gseas <- gseaRes()
    pipelines <- seq(1, length(gseas))
    print(pipelines)
    selectInput("npipeline", "Choose Pipeline", choices=pipelines)
  })
  output$gseaControls <- renderUI({
    req(input$runGsea)
    gseas <- gseaRes()
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
    gseas <- gseaRes()
    print(gseas[[as.numeric(input$npipeline)]])
    print(gseas[[as.numeric(input$npipeline)]][[input$nclusts]])
    plotGeneSetEval(gseas[[as.numeric(input$npipeline)]][[input$nclusts]], clustname=input$nclusts)
  })
  output$metricsControls <- renderUI({
    gseas <- gseaRes()
    pipelines <- seq(1, length(gseas))
    print(pipelines)
    selectInput("npipelineMetrics", "Choose Pipeline", choices=pipelines)
  })
  output$metrics <- renderTable({
    req(input$npipelineMetrics)
    req(input$computeMetrics)
    clusts <- clustRes()
    metrics <- computeUnsupervisedMetrics(sce=readRDS(input$sce$datapath), clusters = clusts[[as.numeric(input$npipelineMetrics)]])
    print(metrics)
    metrics
  }, digits=10,rownames=TRUE, colnames=FALSE)
  
  output$sceSummary <- renderTable({
    req(input$sce)
    sce <- readRDS(input$sce$datapath)
    c(class="SingleCellExperiment", ncells = dim(sce)[[2]], ngenes=dim(embryo)[[1]])
  }, rownames=TRUE, colnames=FALSE)


}
shinyApp(ui, server)