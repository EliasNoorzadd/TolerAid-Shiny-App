# loading in required shiny libraries
library(shiny)
library(shinythemes)
library(shinydashboard)

# the following section defines the user interface of the app
ui <- dashboardPage(
  dashboardHeader(title = "TolerAid"),
  skin = "blue",
  dashboardSidebar(
    
    # creating the menu which is used to switch between pages
    sidebarMenu(
      menuItem("Tolerance Predictor", tabName = "TolerancePredictor", icon = icon("calculator")),
      menuItem("Patient Database", tabName = "PatientDatabase", icon = icon("database"))
    )
  ),
  dashboardBody(
    tabItems(
      # First page
      tabItem(
        tabName = "TolerancePredictor",
        fluidRow(
          tags$h1(style = "margin-left: 15px;", "Tolerance Predictor"),
          # the below box contains the input section for the app
          box(
            title = "Patient Info",
            width = 12,
            # status = "primary",
            color = "black",
            collapsible = TRUE,
            textInput("givenName", label = "Given Name:"),
            textInput("surname", label = "Surname:"),
            dateInput("dob", label = "Date Of Birth:"),
            fileInput("file1", "Choose CSV File",accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
            checkboxInput("immunotherapy","On Immunotherapy"),
            div(
              actionButton("predict","Predict For Sample Patient"),
              style = "position: absolute; bottom: 10px; right: 10px;"
            )
            
          ),
          # here the output of the results and recommendation is outputed
          uiOutput("recommendationBox")
        ),
        # here the output for the boxplot in the more info section is outputed
        fluidRow(
          column(
            width = 12,
            plotOutput("geneBoxplot"),
            textOutput("geneBoxplotExplanation")
          )
        )
      ),
      
      # Second page containing the recent patient database
      tabItem(
        tabName = "PatientDatabase",
        fluidRow(
          #column(
          #tableOutput("showTable"),
          tableOutput("patientDatabase")
          #width = 12
          #)
        )
      )
    )
  )
)

# below the server logic is defined
server <- function(input, output) {
  
  # initialising the contents of the patient table as well as rendering the initial table
  patientTable <- data.frame(givenName = character(), surname = character(), dob = character(), immunotherapyStatus = character(), prediction_result = character(), urgency = character())
  output$patientDatabase <- renderTable(patientTable)
  colnames(patientTable) <- c("Given Name", "Surname", "Date of Birth", "Immunotherapy Status", "Prediction Result", "Urgency")
  
  # loading in required packages for backend logic of the app
  library(tidyverse)
  library(stats)
  library(ggplot2)
  library(randomForest)
  library(here)
  
  currentDirectory <- here()
  currentDirectory
  
  # setting files paths for preprocessed data
  
  boxplotObjectsFilePath <- file.path(currentDirectory, "data/listofobjectsfornewboxplot.RData")
  immunoModelFilePath <- file.path(currentDirectory, "data/rf_model_immunosuppressed.RData")
  notImmunoModelFilePath <- file.path(currentDirectory, "data/rf_model__not_immunosuppressed.RData")
  listOfObjectsFilePath <- file.path(currentDirectory, "data/list_of_objects_for_shiny_app.RData")
  
  load(boxplotObjectsFilePath)
  load(immunoModelFilePath)
  load(notImmunoModelFilePath)
  load(listOfObjectsFilePath)
  
  # the code located within the section below runs once everytime the predict button is clicked
  observeEvent(input$predict,{
    
    # this section of code is responsible for creating the boxplot using the preprocesed gene info
    createnewgeneboxplot=function(immunosuppression_status, patient_name, patient_data){
      if(immunosuppression_status){
        #read in data that goes into model 1 (top 50 genes)
        genedata1=listofobjectsfornewboxplot[[1]]
        top50genelist1=unique(as.character(genedata1$gene))
        
        #make sure to only consider the relevant genes from the input patient
        patient_data_subset1=patient_data[top50genelist1]
        
        #if patient doesnt have all the genes then spit out an error
        if(NA %in% patient_data_subset1){ return(list("Error: Patient sample does not contain the necessary genes. Please ensure that your input csv contains the following genes: ",top50genelist1))}else{
          # find mean expressions of each gene from all the training data samples
          genemeans1=genedata1 %>% group_by(gene) %>%summarise(avg = mean(measurement)) 
          genemeans1$patientvalue=patient_data_subset1
          genemeans1$difference=abs(genemeans1$avg-genemeans1$patientvalue)
          
          #sort list of means by their difference from our sample patient 
          sortedgenemeans1=genemeans1[order(genemeans1$difference, decreasing=TRUE),]
          
          #find the 5 genes with the biggest difference
          biggest_diff_genes=sortedgenemeans1$gene[1:5]
          biggest_diff_data1=genedata1[genedata1$gene %in% biggest_diff_genes,]
          
          #subset patient data to only those 5 genes
          patient_data_biggest_diff1=patient_data[names(patient_data)%in%biggest_diff_genes]
          
          #plot training data
          dataset1plot=ggplot()+geom_boxplot(data=data.frame(biggest_diff_data1), aes(x=measurement, y=gene))
          
          #format subsetted patient data
          patient_df1=data.frame(patient_data_biggest_diff1)
          patient_df1$gene=row.names(patient_df1)
          patient_df1$measurement=patient_df1$patient_data_biggest_diff1
          
          #plot patient data
          exampleplot1=dataset1plot+geom_point(data=patient_df1, aes(x=measurement, y=gene,colour="red",size=1))+ scale_size(guide = 'none')+theme(legend.title=element_blank())+scale_colour_discrete(labels=c(patient_name))+ylab(label="Gene")+xlab(label="Relative Expression")+ggtitle("5 Genes Where Difference Between \n Patient Expression and Average Expression Is Greatest")
          exampleplot1
        }
      }else{
        #exactly the same as the last part but for the non-immunosuppressent dataset
        genedata2=listofobjectsfornewboxplot[[2]]
        top50genelist2=unique(as.character(genedata2$gene))
        patient_data_subset2=patient_data[top50genelist2]
        if(NA %in% patient_data_subset2){ return(list("Error: Patient sample does not contain the necessary genes. Please ensure that your input csv contains the following genes: ",top50genelist2))}else{
          genemeans2=genedata2 %>% group_by(gene) %>%summarise(avg = mean(measurement)) 
          genemeans2$patientvalue=patient_data_subset2
          genemeans2$difference=abs(genemeans2$avg-genemeans2$patientvalue)
          sortedgenemeans2=genemeans2[order(genemeans2$difference, decreasing=TRUE),]
          biggest_diff_genes=sortedgenemeans2$gene[1:5]
          biggest_diff_data2=genedata2[genedata2$gene %in% biggest_diff_genes,]
          patient_data_biggest_diff2=patient_data[names(patient_data)%in%biggest_diff_genes]
          dataset2plot=ggplot()+geom_boxplot(data=data.frame(biggest_diff_data2), aes(x=measurement, y=gene))
          patient_df2=data.frame(patient_data_biggest_diff2)
          patient_df2$gene=row.names(patient_df2)
          patient_df2$measurement=patient_df2$patient_data_biggest_diff2
          
          exampleplot2=dataset2plot+geom_point(data=patient_df2, aes(x=measurement, y=gene,colour="red",size=1))+ scale_size(guide = 'none')+theme(legend.title=element_blank())+scale_colour_discrete(labels=c(patient_name))+ylab(label="Gene")+xlab(label="Relative Expression")+ggtitle("5 Genes Where Difference Between \n Patient Expression and Average Expression Is Greatest")
          exampleplot2
        }
      }
    }
    
    # this reads in the data from the csv input
    file_data <- read.csv(input$file1$datapath)
    geneData <- as.matrix(file_data)[1,]
    output$matrix <- renderDataTable(geneData)
    
    #storing input from patient data
    immunotherapyStatus <- input$immunotherapy
    givenName <- input$givenName
    surname <- input$surname
    dob <- as.character(input$dob)
    predictionProb <- 68
    urgency <- FALSE
    
    # this section is used to make predictions for patients on immunotherapy using the model which is built on patients on immunotherapy as well
    if (immunotherapyStatus == TRUE){
      
      stats::predict(rf_model_immunosuppressed, geneData)
      print("1")
      predictionProbList <- stats::predict(rf_model_immunosuppressed, geneData, type = "prob")
      prediction_result <- as.character(stats::predict(rf_model_immunosuppressed, geneData))
      
      # selecting the confidence outcome for the predicted outcome (aka the one with the largest confidence)
      if (predictionProbList[1]>predictionProbList[2]){
        predictionProb <- predictionProbList[1] * 100
      }
      if (predictionProbList[1]<predictionProbList[2]){
        predictionProb <- predictionProbList[2] * 100
      }
      if (predictionProbList[1] == predictionProbList[2]){
        predictionProb <- 50
      }
      if (prediction_result == "Tolerant") {
        urgency <- TRUE
      }
    }
    
    # This section predicts based on the model built using patients not on immunosuppression
    if (immunotherapyStatus == FALSE){
      stats::predict(rf_model_not_immunosuppressed, geneData)
      predictionProbList <- stats::predict(rf_model_not_immunosuppressed, geneData, type = "prob")
      prediction_result <- as.character(stats::predict(rf_model_not_immunosuppressed, geneData))
      
      # selecting the confidence outcome for the predicted outcome (aka the one with the largest confidence)
      if (predictionProbList[1]>predictionProbList[2]){
        predictionProb <- predictionProbList[1] * 100
      }
      if (predictionProbList[1]<predictionProbList[2]){
        predictionProb <- predictionProbList[2] * 100
      }
      if (predictionProbList[1] == predictionProbList[2]){
        predictionProb <- 50
      }
      if (prediction_result == "Non-Tolerant") {
        urgency <- TRUE
      }
    }
    
    # This section renders the color-coded box for the results section
    output$recommendationBox <- renderUI(
      box(
        title = "Results",
        width = 12,
        
        if (prediction_result == "Tolerant"){
          valueBox(
            "Tolerant", paste0("Confidence: ",as.character(predictionProb),"%"), icon("heart-circle-check", lib = "font-awesome"),
            color = "green", width = 6
          )
        },
        if (prediction_result == "Non-Tolerant"){
          valueBox(
            "Non-Tolerant", paste0("Confidence: ",as.character(predictionProb),"%"), icon("heart-circle-xmark", lib = "font-awesome"),
            color = "red", width = 6
          )
        },
        if (prediction_result == "Stable Recipients On Standard Immunotherapy"){
          valueBox(
            "Stable", paste0("Confidence: ",as.character(predictionProb),"%"), icon("heart-circle-minus", lib = "font-awesome"),
            color = "orange", width = 6
          )
        },
        
        # this box renders the recommendation section of the results along with the associated recommendation for the outcome
        box(
          title = "Recommendation",
          solidHeader = TRUE,
          if (prediction_result == "Tolerant"){
            renderText("The patient's kidney appears to resemble that of normal renal functionality, thus it is recommended they do not partake in immunotherapy whilst being continously monitored.")
          },
          if (prediction_result == "Non-Tolerant"){
            renderText("The patient's kidney appears to resemble that of abnormal renal functionality, thus it is recommended they partake in immunotherapy whilst being continously monitored to prevent kidney failure.")
          },
          if (prediction_result == "Stable Recipients On Standard Immunotherapy"){
            renderText("The patient's kidney appears to resemble that of a stable transplant renal functionality, thus it is recommended they continue with their immunotherapy and should be continously monitored.")
          },
          width = 6
        ),
        
        # this section contains a box for the boxplot to be rendered in the results section
        box(
          title = "More Info",
          width = 12,
          solidHeader = FALSE,
          collapsible = TRUE,
          collapsed = TRUE,
          renderPlot(createnewgeneboxplot(immunotherapyStatus, givenName, geneData)),
          renderText("Above is a boxplot which contains the relative expression of the top 5 genes where the patients expression value significantly varies from the average kidney transplant patient with the same immunotherapy status. This can be used to identify what genes make this specific patient unique from the population.")
        )
      )
    )
    
    # this section renders the table which contains information of the recent patients predicted and their outcomes, as well as the urgency status
    library(stringr)
    cleaned_prediction_result <- str_split_i(prediction_result," ",1)
    patientTable <<- rbind(patientTable, data.frame(givenName, surname, dob, immunotherapyStatus, cleaned_prediction_result, urgency))
    output$patientDatabase <- renderTable(patientTable, colnames(patientTable) <- c("Given Name", "Surname", "Date of Birth", "Immunotherapy Status", "Prediction Result", "Urgency"))
  })
}

# this section is used to run the 2 components of the shiny app (ui and server)
shinyApp(ui = ui, server = server)
