library(shiny)
library(ggplot2)

ui <- fluidPage(
  titlePanel("Enrichment result of two EMI binders"),
  fluidRow(
    column(2,
           fileInput("upload1", NULL, accept = c(".csv", ".tsv"), 
            buttonLabel = "Upload1 .csv")),
    column(2, 
           sliderInput("fcTH1", "Fold change threshold",
                       value = 1.5, min = 1.2, max = 6)),
    column(2, 
           sliderInput("pvalueTH1", "pvalue threshold", 
              value = 0.05, min = 0.01, max = 0.2)),
    column(2,
           fileInput("upload2", NULL, accept = c(".csv", ".tsv"), 
            buttonLabel = "Upload2 .csv")),
    column(2, 
           sliderInput("fcTH2", "Fold change threshold",
                       value = 1.5, min = 1.2, max = 6)),
    column(2, 
           sliderInput("pvalueTH2", "pvalue threshold", 
              value = 0.05, min = 0.01, max = 0.2))
  ),
  fluidRow(
    column(4,
           actionButton("cal", "Calculate"))
  ),
  fluidRow(
    column(6,
           plotOutput("plot1")),
    column(6,
           plotOutput("plot2"))
  ),
  fluidRow(
    column(6,
           tableOutput("enrich1")),
    column(6,
           tableOutput("enrich2"))
  )
)

server <- function(input, output, session) {
  
  data1 <- reactive({
    req(input$upload1)
    
    ext <- tools::file_ext(input$upload1$name)
    switch(ext,
      csv = vroom::vroom(input$upload1$datapath, delim = ","),
      tsv = vroom::vroom(input$upload1$datapath, delim = "\t"),
      validate("Invalid file; Please upload a .csv or .tsv file")
    )
  })
  
  data2 <- reactive({
    req(input$upload2)
    
    ext <- tools::file_ext(input$upload2$name)
    switch(ext,
      csv = vroom::vroom(input$upload2$datapath, delim = ","),
      tsv = vroom::vroom(input$upload2$datapath, delim = "\t"),
      validate("Invalid file; Please upload a .csv or .tsv file")
    )
  })
  
  hit1 <- eventReactive(input$cal, {
    (data1() %>%
       filter(log2fc >= log2(input$fcTH1) & 
               pvalue < input$pvalueTH1 &
               Cys > 0))$proteinID
  })
  
  hit2 <- eventReactive(input$cal, {
    (data2() %>%
       filter(log2fc >= log2(input$fcTH2) & 
               pvalue < input$pvalueTH2 &
               Cys > 0))$proteinID
  })
  
  output$plot1 <- renderPlot({
    plotVolcano(data1(), input$fcTH1, input$pvalueTH1)
  }, width = 400, height = 400, res = 96)
  
  output$plot2 <- renderPlot({
    plotVolcano(data2(), input$fcTH2, input$pvalueTH2)
  }, width = 400, height = 400, res = 96)
  
  output$enrich1 <- renderTable({
    enrich1(hit1())
  })
  
  output$enrich2 <- renderTable({
    enrich1(hit2())
  })
}

shinyApp(ui, server)