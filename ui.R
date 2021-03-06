## ui.R ##
library(DT)
library(shiny)
#library(plotly)
library(shinydashboard)


header <- dashboardHeader(title = "Community Visualization v2.0", titleWidth = 450)

sidebar <- dashboardSidebar(
  p("Communities are ", span("blue", style = "color:#2A9FD6")), 
  p("Entities are ",  span("green", style = "color:#FF8800")),
 # p("Chemicals are ", span("orange", style = "color:#FF8800")), 
  #p(" are ", span("red", style = "color:#CC0000")),
 selectInput("select", label = h5("Select algorithm"), 
             choices = list("Louvain" = "lv", "Walktrap" = "wk", "Fast Greedy" = "fg","Infomap" = "imap","Edge betweeness" = "ebetweens","Label Propagation"="lp","Spinglass"="sg"), 
             selected = "lv"),
 hr(),
  actionButton("back_button", "Back"),
  actionButton("reset_button", "Reset"),  
 radioButtons("interactions", "Show interactions:",
              c("All" = "all",
                "Protein-Protein" = "Protein-Protein",
                "Protein-Disease" = "Protein-Disease",
                "Protein-Chemical" = "Protein-Chemical",
                "Chemical-Disease" = "Chemical-Disease")),
 

 #fluidRow(column(3, verbatimTextOutput("value"))),
 # checkboxGroupInput("node_types", "Entities:",
#                     choices = c("Protein" , "Disease", "Chemical"),
#                     selected = c("Protein" , "Disease", "Chemical")),
  
  textInput("searchentitiy","Search Entity"),
  actionButton("search_button","Search")
)

body <- dashboardBody(
  tags$head(
    tags$script(src='lib/sigma.min.js'),
    tags$script(src='lib/sigma.layout.forceAtlas2.min.js'),
    tags$script(src='lib/sigma.parsers.json.min.js'),
    tags$script(src='rendergraph.js'),
    tags$link(rel = "stylesheet", type = "text/css", href = "graph.css")
  ),
  
  fluidRow(  
    box(     textOutput("name"), 
             uiOutput("graph_with_sigma"),
             title = "Network",
             header = TRUE,
             tags$canvas(id="graph", # graphical output area
                         width="1000",
                         height="800"),tags$div(id="graph2")
    ),
    
    tabBox( title = "Details", 
            id = "details",
            selected = "Entities",
            tabPanel("Entities", DT::dataTableOutput("degree_table")),
            #tabPanel("Degrees", plotlyOutput("degree_distribution")),
            #tabPanel("PageRanks", plotlyOutput("pagerank_distribution")),
            tabPanel("Disease Pathway Info.", plotOutput("pathway_distribution"))
    )
  )
  
)

dashboardPage(header, sidebar, body, skin = "blue")
