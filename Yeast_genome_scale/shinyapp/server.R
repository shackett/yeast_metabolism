setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")


library(shiny) # HTML interface
library(gplots) # textPlot
library(ggplot2) # ggplotting
library(gridExtra) # grid of ggplot plots


# find all pathways
# point each pathway to all relevent reactions using a list - order as factors with ALL being first
# generate figures and dump out when appropriate

##### Five files loaded #####
# pathwaySet - number of members in various pathways
# rxToPW - which pathways is a reaction associated with + name
# reactionInfo - reactiion -> rxMechanism + sanitized names
# pathway_plot_list - plots showing information (performance) at a pathway level
# shiny_flux_data - rxMechanism specific plots

test <- T
if(test == T){
  load("shinyapp/shinySubData.Rdata") # this is a reduced environment for development purposes
}else{
  load("shinyapp/shinyData.Rdata") # this is the full environment
}



shinyServer(function(input, output) {
  
  if(!file.exists("rOCAplots")){
    dir.create("rOCAplots")
  }
  
  # pathway selected
  pathways_selected <- reactive({input$pathway})
  
  # reactions available
  reactions_available <- reactive({sort(rxToPW$reactionName[rxToPW$pathway == pathwaySet$pathway[pathwaySet$display == pathways_selected()]])})
  output$reactions_available <- renderUI({selectInput("reaction_chosen", "Reaction:", as.list(reactions_available()))})
  
  #viewType = reactive({ifelse(input$reaction_chosen == "PATHWAY INFORMATION", "PW", "RX")})
  
  ### If a pathway is chosen ###
  
  PWplots <- reactive({names(pathway_plot_list[[pathways_selected()]])})
  
  output$pw_check <- renderUI({selectInput("pathway_plots", "Pathway plots to choose", choices = as.list(PWplots()))})
  
  chosenPWplots <- reactive({pathway_plot_list[[pathways_selected()]][names(pathway_plot_list[[pathways_selected()]]) %in% input$pathway_plots]})
  
  pw_plot_track <- 0
  output$PW <- renderPlot({
    if(length(chosenPWplots) != 0){
      pwsubList <- chosenPWplots()
      do.call(grid.arrange,  pwsubList)
      
      if(input$pw_save != pw_plot_track){ # incrimenting pw_save indicates that plots should be saved
        pw_plot_track <- input$pw_save
        
        name <- paste0("rOCAplots/", input$pw_filename, ".pdf")
        pdf(file = name, height = 15, width = 15)
        do.call(grid.arrange,  pwsubList)
        dev.off()
        }
      
    }else{
      return(NULL) 
    }
  })
  
  
  ### If a reaction is chosen ###
  
  rID = reactive({rxToPW$rID[rxToPW$reactionName == input$reaction_chosen][1]})
  subtypeChoices <- reactive({reactionInfo$Name[reactionInfo$reaction == rID()]})
  
  output$subforms_available <- renderUI({selectInput("subform_chosen", "Subtype:", as.list(subtypeChoices()))})
  currentRx <- reactive({reactionInfo$rMech[reactionInfo$reaction == rID()][reactionInfo$Name[reactionInfo$reaction == rID()] == input$subform_chosen]})
  
  RXplots <- reactive({names(shiny_flux_data[[currentRx()]]$plotChoices)})
  
  output$rx_check <- renderUI({checkboxGroupInput("reaction_plots", "Reaction plots to choose", choices = as.list(RXplots()))})
  
  chosenRXplots <- reactive({shiny_flux_data[[currentRx()]]$plotChoices[names(shiny_flux_data[[currentRx()]]$plotChoices) %in% input$reaction_plots]})
  
  column_number <- reactive({as.numeric(input$col_num)})
  
  output$RXinfo <- renderPlot({
    textplot(shiny_flux_data[[currentRx()]]$reactionInfo, cex = 3, valign = "top", halign = "left")
    })
  
  #output$test <- reactive({length(chosenRXplots())})
  
  rxn_plot_track <- 0
  output$RX <- renderPlot({
  
    if(length(RXplots) != 0){
      rxsubList <- chosenRXplots()
      rxsubList$ncol = ifelse(length(chosenRXplots()) == 1, 1, column_number())
      
      do.call(grid.arrange,  rxsubList)
      
      if(input$rxn_save != rxn_plot_track){ # incrimenting rxn_save indicates that a plot should be saved
        rxn_plot_track <- input$rxn_save
        
        figDimensions <- as.numeric(ifelse(length(chosenRXplots()) == 1, 15, 25))
        
        name <- paste0("rOCAplots/", input$rxn_filename, ".pdf")
        pdf(file = name, height = figDimensions, width = figDimensions)
        do.call(grid.arrange,  rxsubList)
        dev.off()
       
        }
      
    }else{
      return(NULL) 
    }
  })
  
  
  
})
