library(shiny)
library(ggplot2)

ui <- fluidPage(
  titlePanel("Method Robustness"),
  sidebarLayout(
    sidebarPanel(
      fileInput("upload", NULL, accept = c(".csv", ".tsv"), 
            buttonLabel = "Upload .csv"),
      sliderInput("fcTH", "Fold change threshold", 
              value = 1.5, min = 1.2, max = 6),
      sliderInput("pvalueTH", "pvalue threshold", 
              value = 0.05, min = 0.01, max = 0.2),
      actionButton("cal", "Calculate"),
      plotOutput("plot", brush = "plot_brush"),
      tableOutput("selectedPoints"),
    ),
    mainPanel(
      tableOutput("enrich")
      )
    )
)

server <- function(input, output, session) {
  
  data <- reactive({
    req(input$upload)
    
    ext <- tools::file_ext(input$upload$name)
    switch(ext,
      csv = vroom::vroom(input$upload$datapath, delim = ","),
      tsv = vroom::vroom(input$upload$datapath, delim = "\t"),
      validate("Invalid file; Please upload a .csv or .tsv file")
    )
  })
  
  hit <- eventReactive(input$cal, {
    (data() %>%
      filter(log2fc >= log2(input$fcTH) & 
               pvalue < input$pvalueTH &
               Cys > 0))$proteinID
  })
  
  output$plot <- renderPlot({
    plotVolcano(data(), input$fcTH, input$pvalueTH)
  }, width = 400, height = 400, res = 96)
  
  output$selectedPoints <- renderTable({
    req(input$plot_brush)
    brushedPoints(data(), input$plot_brush)
  })
  
  output$enrich <- renderTable({
    enrich1(hit())
  })
}

shinyApp(ui, server)