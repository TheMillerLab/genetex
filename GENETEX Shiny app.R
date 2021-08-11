##########################################################################################################################
# Load Packages
##########################################################################################################################
library(shiny)
library(tidyverse)
library(shinydashboard)
library(genetex)

##########################################################################################################################
# Create Input Objects
##########################################################################################################################

# Create Platform Dropdown
platform <- c("MGH/SNaPshot","BWH/Oncopanel","Guardant 360","FoundationOne CDx","TEMPUS xT Gene Panel","MSK-IMPACT")

# Create Tissue-type platform
genomics_tissue_type <- c("Primary Cutaneous Tumor",
                          "Metastases",
                          "MCC of Unknown Primary (non-cutaneous lesion at initial presentation)",
                          "Local Recurrence",
                          "Blood/Liquid Biopsy",
                          "Unknown/Not Reported"
                          )

##########################################################################################################################
# Shiny User Interface
##########################################################################################################################
ui <- shinydashboard::dashboardPage(
  ################################################ Title ###################################################
  shinydashboard::dashboardHeader(title = "GENETEX To REDCap"),
  ########################################### Sidebar Content ##############################################
  shinydashboard::dashboardSidebar(
    shinydashboard::sidebarMenu(
      menuItem("Inputs", tabName = "inputs", icon = icon("copy")),
      menuItem("Report", tabName = "report", icon = icon("scroll")),
      menuItem("Data", tabName = "data", icon = icon("th"))
    )
  ),
  ########################################### Page(s) Content ##############################################
  shinydashboard::dashboardBody(
    tabItems(
      #First tab content
      tabItem(tabName = "inputs",
              ###########################################  Tab 1 ##############################################
              fluidRow(
                column(width = 4, offset = 0, style='padding-left:30px; padding-right:30px; padding-top:30px; padding-bottom:30px;',
                       textAreaInput(inputId = "data", label = tagList(icon("paste"),"Genomics Report"), placeholder = "Paste Genomics Report (Required)"),
                       textInput(inputId = "record_id", label = tagList(icon("id-card"), "Record ID"), placeholder = "Enter Subject's Record ID (Required)"),
                       numericInput(inputId = "instrument_instance", label = tagList(icon("window-restore"),"Select Instrument Instance of Report"), value = 1, min = 1)
                       ),
                column(width = 4, offset = 0, style='padding-left:30px; padding-right:30px; padding-top:30px; padding-bottom:30px;',
                       selectInput(inputId = "platform", label = tagList(icon("th-large"),"Select the Genomics Report Platform"), choices = platform),
                       textInput(inputId = "lesion_tag",label = tagList(icon("tag"),"Lesion Tag"), placeholder = "Lesion Tag of the Tissue Sequenced (Optional)"),
                       selectInput(inputId = "genomics_tissue_type", label = tagList(icon("object-group"),"Select Lesion Type"), choices = genomics_tissue_type)
                       ),
                column(width = 4, offset = 0, style='padding-left:30px; padding-right:30px; padding-top:30px; padding-bottom:30px;',
                       textInput(inputId = "date_collected", label = tagList(icon("calendar-day"), 'Date the Tissue was Obtained'), placeholder = "Format: YYYY-MM-DD (Optional)"),
                       textInput(inputId = "redcap_uri", label = tagList(icon("map-marked"), "REDCap Web Addresss"), placeholder = "Enter Web Address of Your Platform (Required)"),
                       passwordInput(inputId = "redcap_api_token", label = tagList(icon("unlock-alt"), "REDCap API Token"), placeholder = "Enter your REDCap API Token (Required)")
                       )
                ),
              fluidRow(
                column(6,
                       actionButton(inputId = "run", label = "Run GENETEX To REDCap", style = 'height:60px; width: 200px; font-size:120%'),
                       align = "center",
                       style = "margin-bottom: 10px;",
                       style = "margin-top: -10px;",
                       style = "margin-left: 300px;")
                )
              ),
      ###########################################  Tab 2 ##############################################
      tabItem(tabName = "report",
              verbatimTextOutput("genomics_report")
              ),
      ###########################################  Tab 3 ##############################################
      tabItem(tabName = "data",
              tableOutput("gtr_table")
              )
      ))
  )

##########################################################################################################################
# Shiny Server Side Commands
##########################################################################################################################

server <- function(input, output, session) {

 ############################  Reactives for the arguments of the function ################################
 # Data Argument
 data.df <- eventReactive(input$run, {
   data <- readr::read_delim(input$data, delim = "\n")
   names(data)[1] <-  "Results"
   data <- data.frame(data)
 })

# Create a reactive for the genetex_to_redcap function
  gtr <- eventReactive(input$run, {
    genetex::genetex_to_redcap(
      data = data.df(),
      record_id = input$record_id,
      instrument_instance = input$instrument_instance,
      lesion_tag = input$lesion_tag,
      platform = input$platform,
      date_collected = input$date_collected,
      genomics_tissue_type = input$genomics_tissue_type,
      redcap_api_token = input$redcap_api_token,
      redcap_uri = input$redcap_uri
    )
  })


  ############################  Text Output of the Genomics Report ################################
  ### This will output in the "Report" tab
  output$genomics_report <- renderText(input$data)


  ############################  Table Output of the Genomics Report ################################
  ### This will output in the "Data: tab
  output$gtr_table <- renderTable(gtr())


}

##########################################################################################################################
# Run the application
##########################################################################################################################
shinyApp(ui, server)
